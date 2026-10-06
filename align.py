from fasta import read_fasta


def smith_waterman(a, b, match=2, mismatch=-1, gap=-2):
    n, m = len(a), len(b)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    best, bi, bj = 0, 0, 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if a[i - 1] == b[j - 1] else mismatch
            H[i][j] = max(0,
                          H[i - 1][j - 1] + s,
                          H[i - 1][j] + gap,
                          H[i][j - 1] + gap)
            if H[i][j] > best:
                best, bi, bj = H[i][j], i, j

    ra, rb = [], []
    i, j = bi, bj
    while i > 0 and j > 0 and H[i][j] > 0:
        s = match if a[i - 1] == b[j - 1] else mismatch
        if H[i][j] == H[i - 1][j - 1] + s:
            ra.append(a[i - 1]); rb.append(b[j - 1]); i -= 1; j -= 1
        elif H[i][j] == H[i - 1][j] + gap:
            ra.append(a[i - 1]); rb.append("-"); i -= 1
        else:
            ra.append("-"); rb.append(b[j - 1]); j -= 1

    return best, "".join(reversed(ra)), "".join(reversed(rb)), (i, bi), (j, bj)


if __name__ == "__main__":
    genome = read_fasta("data/sars.fasta")
    window = genome[950:1100]
    swap = {"A": "T", "T": "A", "C": "G", "G": "C"}

    q = genome[1000:1050]
    tests = {
        "exact": q,
        "1 mismatch": q[:25] + swap[q[25]] + q[26:],
        "1 deletion": q[:25] + q[26:],
    }
    for name, query in tests.items():
        score, qa, wa, (qs, qe), (ws, we) = smith_waterman(query, window)
        print(f"{name}: score={score}, genome position={950 + ws}")
