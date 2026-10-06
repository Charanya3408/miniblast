import random
import time
from fasta import read_fasta
from index import build_index
from search import search, naive_search

random.seed(1)
genome = read_fasta("data/sars.fasta")

t = time.time()
index = build_index(genome)
build = time.time() - t

swap = {"A": "T", "T": "A", "C": "G", "G": "C"}


def mutate(q):
    i = random.randrange(10, len(q) - 10)
    q = q[:i] + swap[q[i]] + q[i + 1:]
    j = random.randrange(10, len(q) - 10)
    return q[:j] + q[j + 1:]


trials = 5
fast_t = naive_t = agree = 0
for _ in range(trials):
    pos = random.randrange(0, len(genome) - 100)
    q = mutate(genome[pos:pos + 60])
    t = time.time()
    f = search(q, genome, index)
    fast_t += time.time() - t
    t = time.time()
    n = naive_search(q, genome)
    naive_t += time.time() - t
    agree += (f[0] == n[0] and f[1] == n[1])

print(f"index build (one-time): {build:.3f}s")
print(f"seed-and-extend: {fast_t / trials * 1000:.1f} ms/query")
print(f"naive full scan: {naive_t / trials * 1000:.0f} ms/query")
print(f"speedup: {naive_t / fast_t:.0f}x")
print(f"same score and position: {agree}/{trials}")
