#!/usr/bin/env python3
"""Reproduce the finite computations in the paper; Python standard library only.

The B and C generators implement different legal-site trees. The separate
literal checker constructs inversion sequences using forbidden triples only.
No generated counts from one method are inputs to another method.
"""
import argparse
from collections import Counter
import csv
import json
from math import comb
from pathlib import Path

from class207_recurrence import counts as recurrence_counts

ROOT = Path(__file__).resolve().parent
SETS = {
    'A152': {'010', '100', '102', '210'},
    'B152': {'011', '201', '210'},
    'B207': {'101', '110', '120', '210'},
    'C207': {'100', '101', '110', '210'},
}


def formula152(n):
    return 1 if n == 0 else (1 + sum(comb(2*j, j) for j in range(n))) // 2


def b_tree(N):
    """Four-type tree; states are (type, slack, maximum gap)."""
    totals = [1]
    states = Counter({('R', 1, 0): 1})
    for n in range(1, N + 1):
        totals.append(sum(states.values()))
        children = Counter()
        for (kind, slack, gap), multiplicity in states.items():
            if kind in ('R', 'U', 'W'):
                children['R', slack+1, 0] += multiplicity
            if kind != 'R':
                stop = gap if kind == 'W' else gap + 1
                for k in range(1, stop):
                    children['V', slack+1, k] += multiplicity
            for j in range(1, slack+1):
                children['W' if kind == 'V' else 'U', slack-j+1, j] += multiplicity
        states = children
    return totals


def c_tree(N):
    """Site-word tree; p records seen S and unused U admissible values."""
    totals = [1]
    states = Counter({('S', 1, 1): 1})
    for n in range(1, N + 1):
        totals.append(sum(states.values()))
        children = Counter()
        for (p, slack, h), multiplicity in states.items():
            for j in range(1, slack+1):
                children[p+'U'*(j-1)+'S', slack-j+1, 1] += multiplicity
            if h:
                children['S', slack+1, 1] += multiplicity
            for i in range(len(p) - h):
                children['U'*p[i+1:].count('U'), slack+1, 0] += multiplicity
        states = children
    return totals


def triple(a, b, c):
    alphabet = sorted({a, b, c})
    return ''.join(str(alphabet.index(v)) for v in (a, b, c))


def legal(word, v, forbidden):
    return all(triple(word[i], word[j], v) not in forbidden
               for j in range(1, len(word)) for i in range(j))


def literal_checks(N=8):
    result = {}
    for name, forbidden in SETS.items():
        words = [()]
        totals = [1]
        checked = 0
        maxima = {}
        slack4 = None
        for n in range(1, N+1):
            words = [w+(v,) for w in words for v in range(n)
                     if legal(w, v, forbidden)]
            totals.append(len(words))
            if name.endswith('152'):
                assert len(words) == formula152(n), (name, n)
                continue
            maxima[n] = Counter(max(w) for w in words)
            if n == 4:
                distribution = Counter(n-max(w) for w in words)
                slack4 = [distribution[s] for s in range(1, 5)]
            for w in words:
                M = max(w)
                values = [v for v in range(n+1) if legal(w, v, forbidden)]
                p = [v for v in values if v <= M]
                old = ''.join('S' if v in w else 'U' for v in p)
                h = int(M in p)
                for v in values:
                    child = w+(v,)
                    childmax = max(child)
                    actual = ''.join('S' if a in child else 'U'
                                     for a in range(childmax+1)
                                     if legal(child, a, forbidden))
                    if v > M:
                        old_kept = old if name == 'C207' else ('S' if h else '')
                        expected = old_kept + 'U'*(v-M-1) + 'S'
                    elif v == M:
                        expected = 'S'
                    else:
                        k = sum(a > v and a not in w for a in p)
                        expected = ('S' if name == 'B207' else '') + 'U'*k
                    assert actual == expected, (name, w, v, actual, expected)
                    expected_slack = n-M-(v-M)+1 if v > M else n-M+1
                    expected_h = int(v >= M)
                    assert n+1-childmax == expected_slack
                    assert int(legal(child, childmax, forbidden)) == expected_h
                    checked += 1
        record = {'counts_0_to_8': totals}
        if name.endswith('207'):
            expected_slack = [6, 9, 5, 1] if name == 'B207' else [6, 10, 4, 1]
            assert slack4 == expected_slack
            record.update(site_transitions_checked=checked, slack_at_length_4=slack4)
        if name == 'C207':
            # Listed R_m numerator coefficients, read in ascending z degree.
            numerators = {
                0: {1: 1}, 1: {2: 1, 3: 1, 4: -1},
                2: {3: 2, 4: 4, 5: -7, 6: 1, 7: 1},
                3: {4: 6, 5: 12, 6: -34, 7: 11, 8: 11, 9: -4, 10: -1},
            }
            for n in range(1, N+1):
                for m, numerator in numerators.items():
                    coefficient = sum(u*comb(n-i+m, m) for i, u in numerator.items() if i <= n)
                    assert maxima[n][m] == coefficient, (n, m)
            record['fixed_maximum_checks'] = {'m': [0, 1, 2, 3], 'lengths': [1, N]}
        result[name] = record
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-data', action='store_true',
                        help='regenerate counts.csv before verifying it')
    args = parser.parse_args()
    B, C, rec = b_tree(12), c_tree(12), recurrence_counts(12)
    assert B == C == rec
    rows = [dict(n=n, A152=formula152(n), B152=formula152(n), B207=B[n], C207=C[n])
            for n in range(13)]
    if args.write_data:
        with (ROOT/'counts.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    with (ROOT/'counts.csv').open(newline='') as stream:
        archived = [{k: int(v) for k, v in row.items()} for row in csv.DictReader(stream)]
    assert archived == rows, 'Archived table mismatch'
    literal = literal_checks()
    for name, result in literal.items():
        assert result['counts_0_to_8'] == [row[name] for row in rows[:9]]
    print(json.dumps({'status': 'PASS', 'tree_and_recurrence_lengths': [0, 12],
                      'literal': literal}, indent=2))


if __name__ == '__main__':
    main()
