# Rule 91: crossings with a bare defect, and conserved quantities

`scat.py` is a solver for the Goursat recurrence in the light-cone
coordinates (module). All programs are in Python 3 without further packages.
All running times are for one core and under 2 seconds.

| statement | program | command | expected output |
|---|---|---|---|
| among the 778 objects with core of length ≤ 12, the bare defect crosses exactly the 41 that are 3-sparse unchanged at both parities (Goursat recurrence) | `t91b.py` (uses `t91a.py`) | `python3 t91b.py 12` | `(1) cores of Phi_91 of length <= 12 : 778  transparent to the bare defect: 41  disagreements with 3-sparse: 0 []` |
| the same, by direct simulation of the cells | `t91c.py` | `python3 t91c.py 12` | `cores of length <= 12 : 778 transparent: 41 disagreements with 3-sparse: 0` |
| the two ray lemmas for 91, on pseudo-random data | `t91b.py` | as above | `(2) first-ray lemma: cases 4000 violations 0 ; second-ray parity lemma: cases 1995 violations 0` |
| crossings that are not free (the object with core 00) | `s91.py` | `python3 s91.py` | counts of the crossings of all pairs of cores of length ≤ 7 with shape kept, shape changed, tail changed, and examples; among them `shape changed ((0,), (0, 0), 0, (4, (0, 0)), (4, (0, 1, 0)), 4, 4)` |
| the densities of Q_1 and Q_2 solve the linear system for r = 2 | `cq_check.py` | `python3 cq_check.py RULE` | `rule 91 Q1 conserved: True  Q2 conserved: True`; for 94 and 95: `Q1 conserved: False  Q2 conserved: True` |
| the identity in the proof of the lemma on Q_1: regrouping, the two cases for s, invariance of every term; Q_1 on all configurations of period ≤ 8 | `q1_identity.py` | `python3 q1_identity.py` | `rule 91 regroup violations 0 description of s violations 0 termwise invariance violations 0 | periodic configs 87380 Q1 violations 0`; for 94 and 95 the last three counts are not 0 |
| all additive conserved quantities with a density on r adjacent cells | `cq.py` | `python3 cq.py 91 2` | `nullity 7`, `rank trivial 4 rank all h 6 nontrivial 2`, then two densities |
| Q_1, Q_2 on all configurations of spatial period ≤ 7 | `q91.py` | `python3 q91.py RULE` | `rule 91 configs 24844 Q1 violations 0 Q2 violations 0`; for 94 and 95 `Q2 violations 0` and `Q1 violations` 15580 and 13296 |

`cq.py RULE r` solves h(F(w)) − h(w_1…w_r) = g(w_1…w_{r+1}) − g(w_0…w_r),
w ∈ A^(r+2), exactly over the rationals, and prints the dimension of the
space of densities h modulo constants and coboundaries (`nontrivial`).
`q91.py` applies one step of F to every configuration of period ≤ 7 (21844
configurations) and to 3000 pseudo-random ones of periods 8..60.
