from index import seed_hits
from align import smith_waterman


def search(query, genome, index, k=11, pad=20, min_votes=3):
    best = None
    for diag, votes in seed_hits(query, index, k).most_common(5):
        if votes < min_votes:
            break
        start = max(0, diag - pad)
        end = min(len(genome), diag + len(query) + pad)
        score, qa, ga, _, (gs, _) = smith_waterman(query, genome[start:end])
        if best is None or score > best[0]:
            best = (score, start + gs, qa, ga)
    return best


def naive_search(query, genome):
    score, qa, ga, _, (gs, _) = smith_waterman(query, genome)
    return score, gs, qa, ga
