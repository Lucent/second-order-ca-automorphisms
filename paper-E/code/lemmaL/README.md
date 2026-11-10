# Lemma L for rules 91, 94 and 95: checks of the lemmas

Notation: f(a,b,c) = a XOR c XOR N(a,b,c); light-cone coordinates q = t − i,
p = t + i. All running times are for one core.

| statement | program | command | expected output |
|---|---|---|---|
| defect climbing (i)–(iii), all binary words of length 22 (hence all shorter ones), both classes | `climb.c` | `gcc -O2 -o climb climb.c; ./climb 22 91; ./climb 22 94` | `rule 91 L=22 steps 8388608  viol(i)=0 viol(iii)=0 viol(ii)=0`, the same for 94 |
| two-step gains of the climbing lemmas, from the Goursat recurrence | `climb_gains.c` | `gcc -O2 -o climb_gains climb_gains.c; for r in 91 94 95; do ./climb_gains $r 20 6; done` | `min 2-step gain 4`, `3`, `4` for 91, 94, 95, and `monotone violations 0` |
| d'Alembert formula | `dalembert.py` | `python3 dalembert.py` | `91 formula violations: (0, 14400)`, and the same line for 94, 95 and 155 |
| separation of finite perturbations, R = 4, all 4^7 fillings of 7 adjacent cells in every valid background | `separ.py` (needs numpy) | `python3 separ.py 91 7 40 4; python3 separ.py 94 7 40 4` | one line per valid background: `separated by T=40: 16384/16384`, `left S_R after entering: 0`, `firings at t>=T/2: 0`; `max sep time` is 14, 16, 14 for 91 (backgrounds 01, 10, 11) and 18, 16 for 94 (backgrounds 00, 11) |
| window closure with base Φ^+ ∪ Φ^-, R ≤ 4, rules 91 and 94 | `../closure/dc.c` | see `../closure/README.md` | `COMPLETE after … passes` |

`climb L RULE` enumerates the 2^L words Z on the positions 0..L−1 and both
classes of the updated positions: (i) the word is unchanged below the lowest
defect and the lowest defect does not move down, (iii) no firing below it,
(ii) after two steps it has moved up by at least 3 (when the word is long
enough to see this). `climb_gains RULE W G` enumerates the 2^W words and both
parities and prints the least gain of the lowest defect after one and after
two column steps, over the words whose lowest defect is below W − G − 1.
`dalembert.py` compares, for 200 pseudo-random diagrams per rule, the value
Y[q,p] with A(q) ⊕ B(p) ⊕ κ(q) ⊕ (parity of the firings counted by the
formula) at the 72 lattice points of a 12 × 12 block: 14400 points per rule.
`separ.py RULE W T R` evolves every filling for T steps and records the
first time at which every window of radius R is a window of a point of
Φ^+ ∪ Φ^-. Running times: under 1 second each, except `separ.py` (a few
seconds).
