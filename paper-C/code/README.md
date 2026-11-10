# Programs for the computations of the paper

This directory contains the programs for the computations specified in the
appendix "Computations" of the paper, and their outputs (`logs/`). The
sections below follow the appendix: the censuses of commuting block maps; the
computation used in proofs (fixed words of all lengths); the
computer-assisted statements; the confirmations of results proved by hand.

Conventions. A cell (x, y) is the integer 2x + y, where y is the current row
and x the previous one, and F(c)_i = (y_i, x_i XOR f(y_{i-1}, y_i, y_{i+1})).
The shift is sigma(c)_i = c_{i+1}. A block map of radius R is a table with
4^(2R+1) entries, indexed by the window read as a number in base 4, leftmost
cell most significant.

Requirements. Python 3 with NumPy. The programs of the section "Censuses"
that use a SAT solver need python-sat (PySAT, with CaDiCaL 1.5.3 and
MiniCard); they are marked "SAT". All other programs run without it:
`centralizer.py` is imported by many of them for the local rule and the block
map utilities, and its import of PySAT is optional. Two programs are in C
(`gcc`). Every program is run from this directory. Unless a command is given,
a program takes no arguments, and its output is the log of the same name in
`logs/`.

Times are for one core; memory is stated when it exceeds 1 GB. The logs of
programs that use random data are reproducible: the seeds are fixed in the
programs. The logs of the SAT programs contain timings, and the order of the
solutions depends on the solver; the numbers and the tables do not.

## Libraries

| file | content |
|---|---|
| `centralizer.py` | the local rule of F_f and of its inverse (`make_rule`), block maps and their composition (`BlockMap`), the tables of sigma^a F^b (`named_maps`), the pair-graph test of injectivity, and the SAT encoding and enumeration of the commuting block maps (Section 3 of the paper) |
| `classify_r4.py` | commutation of a table with F on all 4^(2R+3) windows, with NumPy, for any R |
| `propsearch.py` | the commutation constraints with the equalities merged by union-find; used by `invertible_sat2.py` and `enum54_R3.py` |
| `hw_ca.py`, `hw_common.py` | the local rule as tables and as functions, for the programs `hw*.py` |
| `hw59_census_lib.py` | exact composition of block maps on all windows, for `hw59_census.py` and `hw59_generators.py` |
| `hw_families.py` | families with a hidden involution, for `hw_flipall.py` |

## Censuses (Appendix A.1)

| statement | program and command | output | time |
|---|---|---|---|
| Theorem 4.3: Cent_1(F_54) has 13 elements (5 maps sigma^a F^b, 2 constants, 6 others, none of them injective). Also the 32 maps of radius 1 commuting with F_90 | SAT: `python3 centralizer.py --rules 54,90 --R 1` | `logs/census_54_90_R1.txt` (the 13 tables are listed) | 1 s |
| Theorem 4.3: Cent_2(F_54) exceeds the cap of 2·10^4 solutions | SAT: `python3 centralizer.py --rules 54 --R 2 --max 20000 --noinj` | `logs/census_54_R2_full.txt`, line `20000 solutions (CAPPED)` | 31 min (enumeration 6 min, classification of the 20 000 tables 25 min) |
| Theorem 4.3: the balanced maps of radius 2 commuting with F_54 are 23, all injective; 13 are sigma^a F^b | SAT: `python3 invertible_sat.py --rules 54 --R 2 --native` (MiniCard cardinality constraints); `analyze_inv.py FILE` describes the 10 others when the first command is given `--out FILE` | `logs/invertible_54_R2.txt` | 15 s |
| Section 4.6: the balanced maps of radius 3 commuting with F_54 number more than 2·10^4 | SAT: `python3 invertible_sat2.py --rule 54 --R 3 --max 20001` (writes the tables to `raw_balanced_54_R3.txt`, 330 MB, not included) | `logs/balanced_54_R3.txt`, which ends with `CAPPED` | 19 min, 1.2 GB |
| Section 4.6: the enumeration of the invertible maps of radius 3 commuting with F_54, with injectivity on the configurations of period at most 5 in the clauses and lazily learned witnesses (diamonds of length 13, periodic collisions of period at most 10) shared among the workers. The enumeration is not complete and the paper gives no count | SAT: `python3 enum54_R3.py run --workers 24` then `python3 enum54_R3.py classify` (writes the directory `results/`) | no log is included | days on 24 cores, not finished |
| Proposition 6.16: Cent_2(F_59) has 1158 elements, 133 of them injective | SAT: `python3 centralizer.py --rules 59 --R 1,2 --out census_59_R12.json` | `logs/census_59_R12.txt`; the file `census_59_R12.json` (all tables) is included | 3 min |
| Proposition 6.16: the 133 injective maps are the maps (i)–(iii) of the proposition | `python3 hw59_census.py` (reads `census_59_R12.json`; exact composition on all windows of radius up to 4) | `logs/hw59_census.txt`: `distinct tables found 133`, `matched solutions 133`, `unmatched solutions 0`, `candidates not in census 0` | 2 min |
| 90R at radius 1 and 2: 2^5 and 2^17 commuting maps, by elimination over F_2, with the SAT enumeration of all solutions as a cross-check | `python3 linear_90.py 1,2`; SAT: `python3 centralizer.py --rules 54,90 --R 1` (32 solutions for 90R) and `python3 centralizer.py --rules 150 --R 1,2 --noinj` (16 and 256 solutions for 150R) | `logs/linear_90_R12.txt`, `logs/census_54_90_R1.txt`, `logs/census_150_R12.txt` | see the note below |
| 90R at radius 2: of the 2^17 commuting maps, 98 222 are balanced and 26 are invertible, all affine | `python3 analyze_90.py 2 all` | `logs/analyze_90_R2.txt`, line `enumeration of all 2^17 maps done [5851s]` | 1 h 40 min |
| Theorem 1.6(c), Proposition 7.3: the dimensions of Cent_R(F_90) and Cent_R(F_150) for R ≤ 3 are 5, 17, 65 and 4, 8, 12 | `python3 rank90_packed.py RULE R` for RULE = 90, 150 and R = 1, 2, 3 | `logs/rank_90_R3.txt`, `logs/rank_150_R123.txt` (lines `kernel dim`) | 5 s each |
| Theorem 1.6(c): the affine maps and the two families of maps through the factor maps u, w span a space of dimension 2^(2R)+1, R ≤ 3; every table is re-checked to commute | `python3 span90_R3.py`; `python3 verify_90_structure.py` (R = 1, 2: the span is the whole solution space) | `logs/span_90_R3.txt`, `logs/verify_90_structure.txt` | 6 s; 3 s |
| Remark 8.8: the commuting block maps of radius R ≤ 3 for the affine rules 170, 102, 153, 150 form spaces of dimension 5, 9, 13; 6, 10, 14; 5, 9, 13; 4, 8, 12 | `python3 affine_dim.py RULE 1 2 3` for RULE = 170, 102, 153, 150 (elimination over F_2, independent of `centralizer.py`) | `logs/affine_dim.txt` (R = 1, 2), `logs/affine_dim_R3.txt` | 1 s for R ≤ 2; 25 s for R = 3 |

Notes. (1) The 20 000 solutions at which the enumeration of Cent_2(F_54) is stopped depend on the solver; in the logged run they contain one map sigma^a F^b and the two constants. (2) `linear_90.py 1,2` prints the elimination for R = 1 and 2 within seconds (kernel dimensions 5 and 17) and the classification of the 32 maps of radius 1; it then enumerates and classifies the 2^17 maps of radius 2 without progress output, which takes more than 25 minutes. `logs/linear_90_R12.txt` is its output up to that point; the classification at radius 2 with progress output is `analyze_90.py 2 all`. (3) The SAT enumeration of all solutions for 90R at radius 2, `python3 centralizer.py --rules 90 --R 2 --max 200000 --noinj`, finds 131072 = 2^17 solutions in 12 minutes; `logs/census_90_R2.txt` is its output up to the summary line of the enumeration. The program then classifies the 131 072 tables, which takes hours and was stopped.

## Computation used in proofs (Appendix A.2)

| statement | program | output | time |
|---|---|---|---|
| Proposition 5.2(3), Corollary 1.2: for each of the 256 rules and each k ≤ 40, the words of length k fixed in every context, as walks in the graph on the four pairs; fixed words exist for the 81 rules of 31 classes at length at most 3 and, among the others, only for 133 and 161, at the lengths 4, 7, 8 and 10, …, 40; the search over all words of {0,1,2,3}^k, k ≤ 4, agrees with the reduction to {0,3}^k | `fixed_words.py` | `logs/fixed_words.txt` | 1 s |
| The 88 classes under reflection and the cellwise complement; the 31 classes (81 rules) with a fixed word of length at most 3, with the word (Table 1). Writes `true_classes.json`, which `trains_true.py` reads | `true_classes.py` | `logs/true_classes.txt` | 1 s |
| Table 1: for the listed word of each of the 31 classes, at least two cells can stand between two copies of it (the set C of Theorem 1.1): the construction acts nontrivially | `trains_true.py` (see below) | `logs/trains_true.txt` | 6 s |

## Computer-assisted statements (Appendix A.3)

| statement | program and command | output | time |
|---|---|---|---|
| Sections 6.2 and 6.5: invariant families with a hidden involution, conditions (W1) and (H1), for windows of length at most 6 and closures of at most 300 words, on the 57 classes without a fixed word of length at most 3 | `python3 hw_search_k6.py` (windows of length 6); `python3 hw_search.py K` for K = 4, 5 (all classes, closures of at most 200 words) | `logs/hw_search_k6.txt`, `logs/hw_search_4.txt`, `logs/hw_search_5.txt` | 10 s, 30 s, 3 min |
| Section 6.2: for 59R and p in {6, 12, 30, 60, 120}, every word of X_A and X_B has a context in which F^p changes it | `hw59_period.py` (20 000 random contexts per word and p) | `logs/hw59_period.txt`: every count is below 20000 | 1 min |
| Lemma 6.13: the overlaps of walls of 59R, at distances 2 and 3, on all words of length at most 9 | `hw59_overlap.py`, `hw59_walls.py` | `logs/hw59_overlap.txt`, `logs/hw59_walls.txt` | 1 s |
| Example 7.11: the frames 111 u 111 and their complements: distinct frames share no segment cell; psi_rho preserves the frames, on all 2^31 words of length 31; d + J + sigma^2 J = psi(d) on all 2^21 words of length 21; J(not d) = J(d) on all 2^19 words | `gcc -O2 -o marker_S3 marker_S3.c && ./marker_S3` | `logs/marker_S3.txt` | 11 min |
| After Proposition 6.24: the variant Psi_2 of radius 6 for 219R commutes with F on all 4^15 windows and is an involution | `gcc -O2 -fopenmp -o rot_psi2 rot_psi2.c && ./rot_psi2` (without `-fopenmp` on one core) | `logs/rot_psi2.txt` | 4 min on one core |
| Remark 7.5: Lambda_M is strongly connected for M ≤ 4 for the affine rules 15, 60, 105, 150, 195, 240 and not for 90 and 165; for rule 90 it is not strongly connected for M ≤ 8 | `python3 window_graph.py` (the eight rules, M ≤ 5, and the rank of the linear window system); `python3 window_graph.py 90 8` | `logs/window_graph.txt`, `logs/window_graph_90.txt` | 5 s, 15 s |
| After Theorem 1.4: the orders o_L of F_L for L ≤ 7 and the periods of the content 0^L for L ≤ 10; the cycle lengths of F_L | `walls54.py`, `cycles54.py` | `logs/walls54.txt`, `logs/cycles54.txt` | 4 s, 3 s |
| Table 2: the spread of a changed cell after 300 steps, 40 random configurations of length 801 | `python3 walls_search.py 54 110 30` | `logs/walls_search.txt` | 2 min |

## Confirmations of results proved by hand (Appendix A.4)

None of these is used in a proof.

### Walls

| statement | program | output | time |
|---|---|---|---|
| Theorem 1.3: Psi_m(sigma) on 3000 random cyclic configurations of length about 60 with trains of 2 to 4 cars | `sym54.py` | `logs/sym54.txt`: `commutes with F: True; sigma^-1 inverts: True; homomorphism: True; nontrivial on 476/3000 samples` | 1 s |
| Theorem 1.1: the construction for the word of each of the 31 classes of Table 1, on 1500 random cyclic configurations each | `trains_true.py` (run `true_classes.py` first) | `logs/trains_true.txt`: `construction verified (commutes, acts nontrivially) for 31 rules`, `not verified: []` | 6 s |
| Lemma 4.1 and Section 4: the 32 rules for which 33 is a wall; the maps D^k; the maps Phi_L(F_L) for L = 1, 2, 3 commute with F | `walls54.py`, `walls54b.py` | `logs/walls54.txt`, `logs/walls54b.txt` | 4 s, 2 s |
| Theorem 1.5: the table of tau differs from the identity on 24 of the 16 384 entries, commutes with F on all 4^9 windows, is injective by the pair-graph test and on the configurations of period at most 10, and is not a map sigma^a F^b. Writes `tau_table.json` | `tau.py` | `logs/tau.txt` | 30 s, 4.5 GB |
| Theorem 1.5: tau o tau is the identity on all 4^13 windows of length 13, and tau D differs from D tau on 147 456 of them | `tau_exact.py` | `logs/tau_exact.txt` | 10 s |
| Proposition 5.2: the word 0330 is fixed by F_133 in every context, and D commutes with F on all 4^11 windows; the same for 161 with 3003 | `verify133.py` | `logs/verify133.txt` | 30 s, 1.1 GB |

### Hidden walls: rule 59R (Section 6.2)

| statement | program | output | time |
|---|---|---|---|
| Equation (5), Lemmas 6.11 and 6.12: the table of the walls X_A, X_B and the flip, in all contexts, for F and F^-1 | `hw59_table.py`; a second implementation, `hw59_walls.py` | `logs/hw59_table.txt`, `logs/hw59_walls.txt` | 1 s |
| Lemma 6.13: (H2) and (H4) on all words of length at most 8 | `hw59_overlap.py` | `logs/hw59_overlap.txt`: `bad 0` | 1 s |
| Theorem 6.14: 200 rings of length 2000 seeded with trains of length up to 6, all m = 2, …, 6 permuted at once: F Psi = Psi F for 60 steps, F^-1 Psi = Psi F^-1, the homomorphism law; the radius 2m+2 when all cells outside the window are re-randomized | `hw59_trains.py` | `logs/hw59_trains.txt` | 10 s |
| Proposition 6.16: the eight generators commute with F on all windows of radius 3 | `hw59_generators.py` | `logs/hw59_generators.txt` | 1 s |
| Proposition 6.16: flips preserve the class of every wall, and the 8 classes occur (all words of length 9) | `hw59_classes.py` | `logs/hw59_classes.txt`: `class changes 0` | 1 s |
| Proposition 6.16: the elements of T_R and T_L composed exactly on windows of radius up to 4 | `hw59_census.py` (see "Censuses") | `logs/hw59_census.txt` | 2 min |

### Time-periodic walls: rules 127R, 179R, 223R (Proposition 6.19)

| statement | program | output | time |
|---|---|---|---|
| the cycles of words, in all 16 contexts, for F and F^-1; words half a period apart differ at position 2; overlap offsets | `hwtp_cycles.py`; a second implementation with (H1) and (H2), `hwtp_walls.py` | `logs/hwtp_cycles.txt`, `logs/hwtp_walls.txt` | 1 s |
| the flip of all walls, composed exactly at radius 4, is an involution and passes the pair-graph test | `python3 hw_flipall.py 127,179,223` (uses `hw_families.py`; `python3 hw_families.py 59,27,127,179,223` gives `logs/hw_families.txt`) | `logs/hw_flipall.txt` | 9 s |
| rings of length 1500 with trains at gaps 5, 6 and 8, walls of other classes mixed in: commutation with F and F^-1, homomorphism law | `hwtp_trains.py` | `logs/hwtp_trains.txt`: 0 failures in the nine lines | 10 s |

### Rule 27R (Theorem 6.20)

| statement | program | output | time |
|---|---|---|---|
| (W1) and (H1) for the 12 words in all 16 contexts, for F and F^-1; (H2) on words of length at most 11 | `hw27_walls.py` | `logs/hw27_walls.txt` | 4 s |
| beta and beta + 1 are the only maps with beta o iota = beta + 1 that satisfy (H3) for some epsilon | `hw27_beta.py` | `logs/hw27_beta.txt`: two solutions | 1 s |
| (H4) on words of length at most 11; rings with trains at gaps 3, 6 and 7 | `hw27_trains.py` | `logs/hw27_trains.txt`: `H4 violations 0`, `fails 0` three times | 10 s |

### Rotating walls: rules 219R and 126R (Theorem 6.23, Propositions 6.24 and 6.25)

| statement | program | output | time |
|---|---|---|---|
| the hypotheses on all contexts, in two implementations: bonds are invariant, bonded cells evolve by pi, the maps phi_s and R_s | `rot_walls.py`, `rot_verify.py` | `logs/rot_walls.txt`, `logs/rot_verify.txt` | 2 s |
| D (219R) and E (126R) commute with F and F^-1 on all 4^7 windows and are injective by the pair-graph test | `rot_DE.py`; also item (6) of `rot_verify.py` | `logs/rot_DE.txt` | 6 s |
| 126R satisfies the hypotheses of Theorem 6.6; Corollary 6.9 does not apply to 219R or 126R | `rot_hidden.py` | `logs/rot_hidden.txt` | 1 s |
| Psi_1, Psi_2, Psi_3 on planted rings | `rot_trains.py`; item (5) of `rot_verify.py` | `logs/rot_trains.txt` | 10 s |

### Markers for 90R and 165R (Theorem 7.9, Lemma 7.10)

| statement | program | output | time |
|---|---|---|---|
| the construction of Lemma 7.10 with N = 2 (three words of length 19; the explicit families are checked to lie in U_L for N ≤ 7), Phi and Phi' on planted random rings: commutation with F, the group law, commutation of the two movers; for 165R, the negative control with the frames 111 u 111 alone | `markers.py` | `logs/markers.txt` | 10 s |
| the identities on all small windows, and Phi, Phi' of Example 7.11 on all cyclic configurations of periods 10 and 11 and on planted rings | `marker_cells.py` | `logs/marker_cells.txt` | 17 min, 8 GB |

### Theorem 2.4 and Lemma 2.3

| statement | program | output | time |
|---|---|---|---|
| Theorem 2.4(1), (2) by exhaustion for all 256 rules and n ≤ 4 | `reach_thmA.py` | `logs/reach_thmA.txt`: `violations: []` | 15 s |
| the construction of Lemma 2.3 with random staircase data for the 240 rules that read c and n ≤ 31 | `reach_stair.py` | `logs/reach_stair.txt`: `rules 240 tests 1440 Lemma 2.3 failures 0 Theorem 2.4 failures 0` | 1 s |

## Files written by the programs

`true_classes.json` (by `true_classes.py`) and `tau_table.json` (by `tau.py`)
are included. `search_k4.json`, `search_k5.json`, `search_k6.json` (by
`hw_search.py`, `hw_search_k6.py`), `results_linear90.json` and
`results_analyze90_R2.json` (by `linear_90.py`, `analyze_90.py`),
`raw_balanced_54_R3.txt` and `raw_balanced_54_R3.json` (by `invertible_sat2.py`)
and the directory `results/` (by `enum54_R3.py`) are written again by the
programs and are not included.
