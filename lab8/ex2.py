"""Global DNA sequence alignment using the Needleman-Wunsch algorithm."""

S1 = "ACCGTGAAGCCAATAC"
S2 = "AGCGTGCAGCCAATAC"

MATCH_SCORE = 1
MISMATCH_SCORE = -1
GAP_PENALTY = -2


def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """Return an optimal global alignment, score, and scoring matrix."""
    if not seq1 or not seq2 or any(base not in "ACGT" for base in seq1 + seq2):
        raise ValueError("Both sequences must contain only A, C, G and T.")

    rows, cols = len(seq1) + 1, len(seq2) + 1
    scores = [[0] * cols for _ in range(rows)]

    for i in range(1, rows):
        scores[i][0] = scores[i - 1][0] + gap
    for j in range(1, cols):
        scores[0][j] = scores[0][j - 1] + gap

    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = scores[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            up = scores[i - 1][j] + gap
            left = scores[i][j - 1] + gap
            scores[i][j] = max(diagonal, up, left)

    aligned1, aligned2 = [], []
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            diagonal = scores[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            if scores[i][j] == diagonal:
                aligned1.append(seq1[i - 1])
                aligned2.append(seq2[j - 1])
                i -= 1
                j -= 1
                continue
        if i > 0 and scores[i][j] == scores[i - 1][j] + gap:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    return "".join(reversed(aligned1)), "".join(reversed(aligned2)), scores[-1][-1], scores


def print_matrix(seq1, seq2, scores):
    print("\nScoring matrix (rows: S1, columns: S2):")
    print(f"{'':>5}{'-':>5}" + "".join(f"{base:>5}" for base in seq2))
    for base, row in zip("-" + seq1, scores):
        print(f"{base:>5}" + "".join(f"{value:>5}" for value in row))


def main():
    aligned1, aligned2, score, matrix = needleman_wunsch(
        S1, S2, MATCH_SCORE, MISMATCH_SCORE, GAP_PENALTY
    )
    markers = "".join("|" if a == b else " " if "-" in (a, b) else "." for a, b in zip(aligned1, aligned2))
    matches = sum(a == b for a, b in zip(aligned1, aligned2))
    mismatches = sum(a != b and "-" not in (a, b) for a, b in zip(aligned1, aligned2))
    gaps = sum("-" in (a, b) for a, b in zip(aligned1, aligned2))

    print("GLOBAL DNA ALIGNMENT — Needleman-Wunsch")
    print(f"S1: {S1}\nS2: {S2}")
    print(f"Scoring: match={MATCH_SCORE}, mismatch={MISMATCH_SCORE}, gap={GAP_PENALTY}")
    print_matrix(S1, S2, matrix)
    print("\nOptimal global alignment:")
    print("S1: " + aligned1)
    print("    " + markers)
    print("S2: " + aligned2)
    print(f"\nAlignment score: {score}")
    print(f"Matches: {matches} | Mismatches: {mismatches} | Gaps: {gaps}")
    print(f"Identity: {100 * matches / len(aligned1):.2f}%")


if __name__ == "__main__":
    main()
