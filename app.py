import hashlib
import re
import time
from collections import OrderedDict
from flask import Flask, render_template_string, request
from markupsafe import Markup
from fasta import read_fasta
from index import build_index
from search import search

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 1_000_000

DEFAULT_ID = "sars"
DEFAULT_NAME = "SARS-CoV-2 (NC_045512.2)"
GENOME = read_fasta("data/sars.fasta")
INDEX = build_index(GENOME)
COMP = str.maketrans("ACGT", "TGCA")
MAX_LEN = 500
MAX_GENOME = 500_000
CACHE = OrderedDict()


def clean(raw):
    lines = [l for l in raw.splitlines() if not l.startswith(">")]
    return re.sub(r"\s+", "", "".join(lines)).upper()


def revcomp(s):
    return s.translate(COMP)[::-1]


def span(s):
    cls = {"A": "a", "C": "c", "G": "g", "T": "t"}
    return Markup("".join(f'<b class="{cls.get(c, "x")}">{c}</b>' for c in s))


def render_alignment(qa, ga, width=60):
    blocks = []
    for i in range(0, len(qa), width):
        q, g = qa[i:i + width], ga[i:i + width]
        mid = "".join("|" if x == y else " " for x, y in zip(q, g))
        blocks.append({"q": span(q), "mid": mid, "g": span(g)})
    return blocks


def make_samples(genome, start):
    s = ""
    for p in range(start, len(genome) - 60, 60):
        w = genome[p:p + 60]
        if not set(w) - set("ACGT"):
            s = w
            break
    if not s:
        return "", ""
    s = s[:20] + ("A" if s[20] != "A" else "T") + s[21:40] + s[41:]
    return s, revcomp(s)


def parse_upload(data):
    text = data.decode("utf-8", errors="ignore")
    seq, name, seen = [], "uploaded genome", False
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if seen:
                break
            seen = True
            name = line[1:].strip()[:60] or name
            continue
        seq.append(line)
    genome = re.sub(r"\s+", "", "".join(seq)).upper()
    genome = re.sub(r"[^ACGT]", "N", genome)
    if len(genome) < 100:
        raise ValueError("Genome is too short (minimum 100 bases) or not a FASTA file.")
    if len(genome) > MAX_GENOME:
        raise ValueError(f"Genome is too long (maximum {MAX_GENOME:,} bases).")
    if genome.count("N") / len(genome) > 0.10:
        raise ValueError("More than 10% of the bases are not A, C, G or T.")
    return name, genome


def add_genome(name, genome):
    gid = hashlib.sha1(genome.encode()).hexdigest()[:12]
    if gid in CACHE:
        CACHE.move_to_end(gid)
    else:
        CACHE[gid] = (name, genome, build_index(genome))
        while len(CACHE) > 1:
            CACHE.popitem(last=False)
    return gid


def run(query, genome, index):
    t = time.time()
    hits = []
    fwd = search(query, genome, index)
    rev = search(revcomp(query), genome, index)
    ms = (time.time() - t) * 1000
    if fwd:
        hits.append((fwd, "forward"))
    if rev:
        hits.append((rev, "reverse complement"))
    if not hits:
        return None, ms
    (score, pos, qa, ga), strand = max(hits, key=lambda h: h[0][0])
    matches = sum(1 for x, y in zip(qa, ga) if x == y)
    span_len = sum(1 for c in ga if c != "-")
    return {
        "score": score, "strand": strand,
        "start": pos + 1, "end": pos + span_len,
        "identity": round(100 * matches / len(qa), 1),
        "length": len(qa), "qlen": len(query),
        "blocks": render_alignment(qa, ga),
    }, ms


@app.errorhandler(413)
def too_large(e):
    return "File too large. Maximum is 1 MB (500,000 bases).", 413


@app.route("/", methods=["GET", "POST"])
def home():
    posted = request.method == "POST"
    raw = request.form.get("query", "")
    gid = request.form.get("gid", DEFAULT_ID)
    f = request.files.get("genome")
    uploaded = posted and bool(f and f.filename)
    result = ms = error = notice = None
    searched = False

    if uploaded:
        try:
            gname, g = parse_upload(f.read())
            gid = add_genome(gname, g)
        except ValueError as e:
            error = str(e)
            gid = DEFAULT_ID

    name, genome, index = DEFAULT_NAME, GENOME, INDEX
    if gid != DEFAULT_ID:
        if gid in CACHE:
            name, genome, index = CACHE[gid]
        else:
            gid = DEFAULT_ID
            error = error or "Uploaded genome expired. Please upload it again."

    if posted and not error:
        q = clean(raw)
        if uploaded and not q:
            notice = "Genome loaded. Now paste a query and press Search."
        elif len(q) < 15:
            error = "Enter at least 15 bases."
        elif len(q) > MAX_LEN:
            error = f"Maximum {MAX_LEN} bases."
        elif set(q) - set("ACGT"):
            error = "Only A, C, G, T are allowed."
        else:
            result, ms = run(q, genome, index)
            searched = True

    start = 5000 if gid == DEFAULT_ID else len(genome) // 2
    sample, sample_rc = make_samples(genome, start)
    return render_template_string(
        PAGE, raw=raw, result=result, ms=ms, error=error, notice=notice,
        searched=searched, glen=len(genome), gname=name, gid=gid,
        default_id=DEFAULT_ID, sample=sample, sample_rc=sample_rc)


PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MiniBLAST</title>
<style>
*{box-sizing:border-box}
body{margin:0;font-family:system-ui,sans-serif;background:#0f172a;color:#e2e8f0}
header{padding:40px 20px 56px;text-align:center;background:linear-gradient(135deg,#6366f1,#ec4899)}
header h1{margin:0;font-size:2.4rem;color:#fff}
header p{margin:8px 0 0;color:#fce7f3}
main{max-width:900px;margin:-28px auto 40px;padding:0 16px}
.card{background:#1e293b;border-radius:14px;padding:20px;margin-bottom:16px;box-shadow:0 8px 24px #0006}
textarea{width:100%;height:120px;background:#0f172a;color:#e2e8f0;border:1px solid #475569;border-radius:8px;padding:10px;font:14px ui-monospace,Menlo,monospace}
button{border:0;border-radius:8px;padding:10px 18px;font-weight:600;cursor:pointer;background:#6366f1;color:#fff}
button.alt{background:#334155;margin-left:6px}
.note{color:#94a3b8;font-size:.92rem;margin:0 0 10px}
.note b{color:#e2e8f0}
a{color:#a5b4fc}
input[type=file]{color:#94a3b8;font-size:.9rem}
.stats{display:flex;gap:12px;flex-wrap:wrap}
.stat{flex:1;min-width:130px;background:#0f172a;border-radius:10px;padding:12px;text-align:center}
.stat b{display:block;font-size:1.3rem;color:#a5b4fc}
.aln{overflow-x:auto;font:14px/1.5 ui-monospace,Menlo,monospace;white-space:pre}
.aln div{margin-bottom:14px}
.aln b{font-weight:600}
.a{color:#4ade80}.c{color:#60a5fa}.g{color:#fb923c}.t{color:#f87171}.x{color:#94a3b8}
.err{color:#fca5a5}
</style></head><body>
<header><h1>MiniBLAST</h1>
<p>Seed-and-extend DNA sequence search</p></header>
<main>
<div class="card"><form method="post" enctype="multipart/form-data">
<input type="hidden" name="gid" value="{{ gid }}">
<p class="note">Genome: <b>{{ gname }}</b> ({{ "{:,}".format(glen) }} bp){% if gid != default_id %} &middot; <a href="/">use SARS-CoV-2 instead</a>{% endif %}</p>
<p class="note">Search your own genome: <input type="file" name="genome" accept=".fasta,.fa,.fna,.txt"><br>FASTA, first record only, max 500,000 bases.</p>
<textarea name="query" id="q" placeholder="Paste a DNA sequence (A, C, G, T), 15 to 500 bases. A FASTA header line is fine.">{{ raw }}</textarea>
<p><button>Search</button>
{% if sample %}<button type="button" class="alt" onclick='fill({{ sample|tojson }})'>Example</button>
<button type="button" class="alt" onclick='fill({{ sample_rc|tojson }})'>Reverse-strand example</button>{% endif %}</p>
</form></div>
{% if error %}<div class="card err">{{ error }}</div>{% endif %}
{% if notice %}<div class="card">{{ notice }}</div>{% endif %}
{% if searched and not result %}<div class="card">No significant match found.</div>{% endif %}
{% if result %}
<div class="card"><div class="stats">
<div class="stat"><b>{{ result.score }}</b>score</div>
<div class="stat"><b>{{ result.identity }}%</b>identity</div>
<div class="stat"><b>{{ "{:,}".format(result.start) }}-{{ "{:,}".format(result.end) }}</b>genome position</div>
<div class="stat"><b>{{ "%.1f"|format(ms) }} ms</b>search time</div></div>
<p>Strand: {{ result.strand }}. Aligned {{ result.length }} columns of a {{ result.qlen }}-base query.</p></div>
<div class="card aln">{% for b in result.blocks %}<div>Query   {{ b.q }}
        {{ b.mid }}
Genome  {{ b.g }}</div>{% endfor %}</div>
{% endif %}
</main>
<script>function fill(s){document.getElementById('q').value=s}</script>
</body></html>"""


if __name__ == "__main__":
    app.run(port=5001, debug=True)
