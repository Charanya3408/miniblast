# MiniBLAST

A seed-and-extend DNA sequence search tool in pure Python, built from scratch. It uses the same core idea as BLAST: index short k-mers, find candidate locations from exact seed matches, then run Smith-Waterman local alignment only around those locations.

## Result

Benchmarked on the SARS-CoV-2 reference genome (NC_045512.2, 29,903 bp) with 30 random 60 bp queries, each carrying one substitution and one deletion:

| Method | Time per query |
|---|---|
| Naive full-genome Smith-Waterman | 203 ms |
| MiniBLAST (seed-and-extend) | 1.4 ms |
| **Speedup** | **149x** |

Both methods returned the same best score and genome position on **30/30** queries. The k-mer index is built once, in about 5 ms.

## How it works

1. **Index** (`index.py`): map every 11-mer in the genome to its positions.
2. **Seed** (`index.py`): look up each 11-mer of the query. Each hit votes for a diagonal (genome position minus query position). A true match piles many votes on one diagonal.
3. **Extend** (`align.py`, `search.py`): run Smith-Waterman on a small window around the top diagonals only, with scoring match +2, mismatch -1, gap -2.

## Run it

```bash
python3 align.py       # alignment tests: exact, 1 mismatch, 1 deletion
python3 benchmark.py   # speed and correctness vs naive scan
```

Needs Python 3 and no external packages. The genome is in `data/sars.fasta`.

## Files

- `fasta.py`: FASTA reader
- `index.py`: k-mer index and seed voting
- `align.py`: Smith-Waterman with traceback
- `search.py`: seed-and-extend search, plus the naive baseline
- `benchmark.py`: timing and agreement check

## Limitations

- The naive baseline is a full Smith-Waterman scan, not the real BLAST program. The speedup is against that baseline only.
- Queries in the benchmark are taken from the genome itself and lightly mutated, so they are easy to find. Divergent or random queries were not tested.
- The search is a heuristic. A query with fewer than 3 seed votes returns no result, which can miss a very divergent match.
- Pure Python, tested on one 30 kb genome only. Larger genomes would need a faster implementation.
