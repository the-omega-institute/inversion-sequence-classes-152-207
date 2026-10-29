#!/usr/bin/env python3
"""Class 207(2): exact descendant recurrence, using integer partial sums."""
import argparse
import csv


def counts(N):
    """Return the counts at lengths 0,...,N; O(N^3) integer operations."""
    if N < 0:
        raise ValueError("N must be nonnegative")
    if N == 0:
        return [1]
    # A[s][k] = a_r(s,k), B[s] = b_r(s), on s+k <= N-r.
    A = [[]] + [[1] * (N - s + 1) for s in range(1, N + 1)]
    B = [0] + [1] * N
    answer = [1, 1]
    for L in range(N, 1, -1):
        # Row sums, diagonal sums, and sums of b_r(t)-a_r(t,0).
        row = [[]]
        for s in range(1, L + 1):
            p = [0]
            for v in A[s]:
                p.append(p[-1] + v)
            row.append(p)
        diag = [[0] for _ in range(L + 1)]
        for m in range(1, L + 1):
            for s in range(1, m + 1):
                diag[m].append(diag[m][-1] + A[s][m - s])
        correction = [0]
        for s in range(1, L + 1):
            correction.append(correction[-1] + B[s] - A[s][0])
        A = [[]] + [
            [row[s + 1][k] + diag[s + k][s] + correction[s]
             for k in range(L - s)] for s in range(1, L)]
        B = [0] + [B[s + 1] + diag[s + 1][s] + correction[s]
                   for s in range(1, L)]
        answer.append(B[1])
    return answer


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("N", nargs="?", type=int, default=12)
    parser.add_argument("--check", metavar="CSV", help="compare n,B207,C207 CSV")
    args = parser.parse_args()
    values = counts(args.N)
    if args.check:
        with open(args.check, newline="") as stream:
            rows = {int(r["n"]): r for r in csv.DictReader(stream)}
        for n, value in enumerate(values):
            if n not in rows or any(int(rows[n][c]) != value
                                    for c in ("B207", "C207")):
                raise SystemExit(f"Reference mismatch at n={n}")
        print(f"Verified all {len(values)} terms (0 <= n <= {args.N}).")
    else:
        print("n,count_207_2")
        for n, value in enumerate(values):
            print(f"{n},{value}")
