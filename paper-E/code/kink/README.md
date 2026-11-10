# Kink rigidity for rules 94 and 95: checks of the lemmas

Cells are c_i = (x_i, y_i), and F(c)_i = (y_i, x_i XOR f(y_{i-1}, y_i, y_{i+1})).
All programs use the direct simulation of F in the coordinates (t, i) of
`sim.py` (needs numpy): the configuration is a point P of Φ^+ on the cells
i < 0 and a point X of Φ^- on the cells i ≥ 0, given by their diagrams π(q),
q = t − i, and ξ(p), p = t + i. `sim.run` returns the row at time T and the
list of firings [q, p]; `sim.outputs` reads off the outgoing parts γ(p) and
o(q). No program uses the Goursat recurrence. All running times are for one
core.

| statement | program | command | expected output |
|---|---|---|---|
| crossing of two kinks for 94, all eight cases | `t_pair.py` | `python3 t_pair.py` | 8 lines (m, a, b, firings, kinks of Γ and Ω): no firing for m = 0 or a ≢ b; for m = 1, a ≡ b the firings (a−1, b−1), (a+1, b+1) and both kinks delayed by 1 |
| one exact ray (94), both front values and the version with q and p exchanged | `t_R1.py` | `python3 t_R1.py 10` | 4 lines `front F mirror M cases N count mismatches 0 location mismatches 0`; N = 266000, 266000, 466930, 466930 (1,465,860 in all); log `logs/t_R1_10.log` |
| d'Alembert formula in the form used in Step 3 | `t_dal.py` | `python3 t_dal.py` | `cases 800 bad 0` |
| train lemma, 2 ≤ g ≤ 8 and g ∈ {15, 22, 31}, at the four placements a, b ∈ {8, 9} (the Φ^+ kink emerges at a, or at a + 1 when a ≡ b + g) | `t_train.py` | `python3 t_train.py` | the lines `beta 0 …` (the background of the lemma) are printed for g = 3 and g = 31 only and end in `ok True`; a failing case would be printed with `ok False`. The lines `beta 1 …` concern the train on the background 1, where `ok False` occurs for g ≤ 2 only |
| Steps 2–5 end to end: 351 words w | `t_step4.py` | `python3 t_step4.py` | `cores 351`, `D0 cores passing both shifts: [()] 1` |
| crossing of two defects for 95 | `t95a.py` | `python3 t95a.py` | 4 lines; two firings and both zeros delayed by 2 for a ≡ b, none for a ≢ b |
| rows for 95 (both parts) and the relations (T), (Δ) on all pairs of cores of length ≤ 9 | `t95b.py` | `python3 t95b.py 9` | `cores 40`, `first-row tests 6000 mismatch 0 second-row tests 6000 mismatch 0`, `solutions of both relations: [((0,), (0,), 11, 10)]` |

`t_R1.py NMAX` enumerates all continuations of length ≤ NMAX of the two
sequences, at two values of a_1 and two lengths of the front of the object.
`python3 t_R1.py 8` (about 25 seconds) prints the case numbers 28080, 28080,
49192, 49192; `python3 t_R1.py 10` takes about 5 minutes. The other programs
take under 2 seconds each. `t95b.py K` tests the two parts of the row lemma
on all pairs of cores of length ≤ K (words of Φ_95 that begin and end with
0) at two placements and both classes (6000 tests each for K = 9).
