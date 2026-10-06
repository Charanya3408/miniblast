def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    return "".join(seq)


if __name__ == "__main__":
    genome = read_fasta("data/sars.fasta")
    print(len(genome))
    print(genome[:50])
