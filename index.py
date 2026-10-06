from collections import Counter, defaultdict
from fasta import read_fasta


def build_index(genome, k=11):
    index = defaultdict(list)
    for i in range(len(genome) - k + 1):
        index[genome[i:i + k]].append(i)
    return index


def seed_hits(query, index, k=11):
    # diagonal = genome position minus query position.
    # A true match produces many hits on the same diagonal.
    diagonals = Counter()
    for q in range(len(query) - k + 1):
        for g in index.get(query[q:q + k], []):
            diagonals[g - q] += 1
    return diagonals


if __name__ == "__main__":
    genome = read_fasta("data/sars.fasta")
    index = build_index(genome)
    print("distinct k-mers:", len(index))

    query = genome[1000:1050]
    print(seed_hits(query, index).most_common(3))
