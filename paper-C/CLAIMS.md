# Paper C: statements and the computations behind them

"Every finite group from a wall: large automorphism groups of second-order
cellular automata". Statement numbers are those of the compiled PDF
(`tex/main.pdf`); the appendix is Appendix A, "Computations". Every command
is run from `code/`; logs are in `code/logs/`. Status:

- **computer-assisted**: the statement, or the part named, rests on the computation;
- **confirmation only**: the statement is proved by hand and the computation checks it; it is not used in a proof.

Tier **quick** is run by `make check` (compared with the log byte for byte,
or, for the SAT programs, on the lines that do not contain timings); tier
**full** is not, and the log listed is the record of the run. Times are for
one core, on the machine described in `environment/README.md`; memory is
below 1 GB unless stated. "SAT" marks the programs that need python-sat.

## Censuses (Appendix A.1)

| # | statement | status | program | command | expected output / log | time, memory | tier |
|---|---|---|---|---|---|---|---|
| 1 | Theorem 4.3, radius 1: Cent_1(F_54) has 13 elements (5 maps σ^a F^b, 2 constants, 6 non-injective others) | computer-assisted | `centralizer.py` (SAT) | `python3 centralizer.py --rules 54,90 --R 1` | `logs/census_54_90_R1.txt`: `13 solutions`, `'named': 5, 'const': 2, 'other': 6, 'invertible': 5`, and the 13 tables | 1 s | quick |
| 2 | Theorem 4.3, radius 2: Cent_2(F_54) exceeds the cap of 2·10^4 solutions | computer-assisted | `centralizer.py` (SAT) | `python3 centralizer.py --rules 54 --R 2 --max 20000 --noinj` | `logs/census_54_R2_full.txt`: `20000 solutions (CAPPED)` | 31 min | full |
| 3 | Theorem 4.3, radius 2: exactly 23 balanced commuting maps, all injective by the pair graph; 13 are σ^a F^b, 10 others | computer-assisted | `invertible_sat.py` (SAT, MiniCard) | `python3 invertible_sat.py --rules 54 --R 2 --native` | `logs/invertible_54_R2.txt`: `balanced commuting maps: 23`, `invertible: 23; of which sigma^a F^b: 13 (expected 13); OTHER invertible: 10`, and the tables | 15 s | quick |
| 4 | Theorem 4.3, radius 2: each of the 10 other maps differs from some F^b on the three windows 33s33 | computer-assisted | `analyze_inv.py` | `python3 invertible_sat.py --rules 54 --R 2 --native --out inv.json && python3 analyze_inv.py inv.json` | printed description of the 10 maps; no log included | 20 s | full |
| 5 | Section 4.6 and Appendix A.1: the balanced commuting maps of radius 3 for 54R number more than 2·10^4 | computer-assisted | `invertible_sat2.py` (SAT, MiniCard) | `python3 invertible_sat2.py --rule 54 --R 3 --max 20001` | `logs/balanced_54_R3.txt`: `20000 balanced solutions`, then `CAPPED` | 19 min, 1.2 GB | full |
| 6 | Section 4.6: the enumeration of Cent^×_3(F_54) with injectivity on periodic configurations and lazily learned witnesses. The enumeration is not complete; the paper uses only qualitative statements from it and gives no count | computer-assisted (incomplete run) | `enum54_R3.py` (SAT) | `python3 enum54_R3.py run --workers 24`, then `python3 enum54_R3.py classify` | writes `results/`; no log is included in this release | not finished after days on 24 cores | full |
| 7 | Proposition 6.16: Cent_2(F_59) was enumerated; 133 of its 1158 elements pass the pair-graph test | computer-assisted | `centralizer.py` (SAT) | `python3 centralizer.py --rules 59 --R 1,2 --out census_59_R12.json` | `logs/census_59_R12.txt`: `1158 solutions`, `'invertible': 133`; `census_59_R12.json` | 3 min | full |
| 8 | Proposition 6.16: the 133 invertible maps are the maps (i)–(iii), built by exact composition on windows of radius up to 4 | computer-assisted | `hw59_census.py`, `hw59_census_lib.py` | `python3 hw59_census.py` | `logs/hw59_census.txt`: `distinct tables found 133 multi: []`, `matched solutions 133`, `unmatched solutions 0 []`, `candidates not in census 0 []` | 2 min | full |
| 9 | Appendix A.1, 90R at radius 2: of the 2^17 commuting maps 98 222 are balanced and 26 invertible, all affine (5851 s) | computer-assisted | `analyze_90.py`, `linear_90.py` | `python3 analyze_90.py 2 all` | `logs/analyze_90_R2.txt`: `enumeration of all 2^17 maps done [5851s]: {'const': 2, 'balanced': 98222, 'other': 131057, 'invertible': 26, 'invertible_affine': 26, 'named': 13}` | 1 h 40 min | full |
| 10 | Theorem 1.6(c): dim Cent_R(F_90) = 5, 17, 65 for R = 1, 2, 3, by elimination | computer-assisted | `rank90_packed.py` | `python3 rank90_packed.py 90 R`, R = 1, 2, 3 | `kernel dim 5`, `17`, `65`; `logs/rank_90_R3.txt` for R = 3 | 5 s | quick (R = 3) |
| 11 | Proposition 7.3: dim Cent_R(F_150) = 4, 8, 12 for R = 1, 2, 3, by elimination | computer-assisted | `rank90_packed.py` | `python3 rank90_packed.py 150 R`, R = 1, 2, 3 | `logs/rank_150_R123.txt`: `kernel dim 4`, `8`, `12` | 5 s | quick (R = 3) |
| 12 | Appendix A.1: for R ≤ 2 the dimensions are cross-checked against the SAT enumeration of all solutions | computer-assisted | `centralizer.py` (SAT), `linear_90.py` | `python3 centralizer.py --rules 54,90 --R 1` (90R: 32 solutions); `python3 centralizer.py --rules 150 --R 1,2 --noinj` (16 and 256 solutions); `python3 linear_90.py 1,2` (elimination, and classification of all 2^5 and 2^17 solutions); `python3 centralizer.py --rules 90 --R 2 --max 200000 --noinj` | `logs/census_54_90_R1.txt`, `logs/census_150_R12.txt`, `logs/linear_90_R12.txt`; `logs/census_90_R2.txt`: `131072 solutions` (the log ends at this line; see `code/README.md`, note 3) | 1 s; 30 s; seconds for the elimination (see `code/README.md`, note 2); 12 min | quick (90R, R = 1), full (the others) |
| 13 | Theorem 1.6(c): the three families span a space of dimension 2^(2R)+1 for R ≤ 3, each table re-verified to commute; for R ≤ 2 the span is the whole solution space | computer-assisted | `span90_R3.py`, `verify_90_structure.py` | `python3 span90_R3.py`; `python3 verify_90_structure.py` | `logs/span_90_R3.txt`: `rank(A+U+W)=5`, `17`, `65`; `logs/verify_90_structure.txt`: `EXACT SPAN: True` twice | 6 s; 3 s | quick |
| 14 | Remark 8.8: for f in {170, 102, 153, 150} and R ≤ 3 the solution spaces have dimensions 5, 9, 13; 6, 10, 14; 5, 9, 13; 4, 8, 12 | computer-assisted | `affine_dim.py` | `python3 affine_dim.py RULE 1 2 3` for RULE = 170, 102, 153, 150 | `logs/affine_dim.txt` (R = 1, 2), `logs/affine_dim_R3.txt` (R = 3) | 1 s for R ≤ 2; 25 s for R = 3 | quick (R ≤ 2), full (R = 3) |

## Computations used in proofs (Appendix A.2)

| # | statement | status | program | command | expected output / log | time, memory | tier |
|---|---|---|---|---|---|---|---|
| 15 | Proposition 5.2(3), Corollary 1.2: fixed words of every length k ≤ 40 for the 256 rules; they exist for the 81 rules of Table 1 at length at most 3 and otherwise only for 133 and 161, at the lengths 4, 7, 8 and 10 ≤ k ≤ 40; the search over all words of length k ≤ 4 agrees with the reduction to {0,3}^k | computer-assisted | `fixed_words.py` | `python3 fixed_words.py` | `logs/fixed_words.txt` | 1 s | quick |
| 16 | Table 1: the 31 classes (81 rules) with a fixed word of length at most 3 | computer-assisted | `true_classes.py` | `python3 true_classes.py` | `logs/true_classes.txt`: `88 true classes; walled rules 81; walled classes 31; wall-free classes 57` | 1 s | quick |
| 17 | Appendix A.2: for each rule of Table 1 the set C of Theorem 1.1 for the listed word has at least two elements | computer-assisted | `trains_true.py` | `python3 true_classes.py && python3 trains_true.py` | `logs/trains_true.txt`: the construction acts nontrivially for the 31 class representatives, `not verified: []` | 6 s | quick |

## Computer-assisted statements (Appendix A.3)

| # | statement | status | program | command | expected output / log | time, memory | tier |
|---|---|---|---|---|---|---|---|
| 18 | Sections 6.2, 6.5: the search for invariant families with a hidden involution satisfying (W1) and (H1), windows of length at most 6, closures of at most 300 words | computer-assisted | `hw_search.py`, `hw_search_k6.py` | `python3 hw_search.py 4`; `python3 hw_search.py 5`; `python3 hw_search_k6.py` | `logs/hw_search_4.txt`, `logs/hw_search_5.txt`, `logs/hw_search_k6.txt` | 10 s, 30 s, 3 min | full |
| 19 | Section 6.2: for 59R and p in {6, 12, 30, 60, 120}, every word of X_A ∪ X_B has a context in which F^p changes the protector | computer-assisted (random search, fixed seed) | `hw59_period.py` | `python3 hw59_period.py` | `logs/hw59_period.txt`: all 24 × 5 counts are below 20000 | 1 min | full |
| 20 | Lemma 6.13: the overlaps of walls of 59R, at distances 2 and 3, on all words of length at most 9 | computer-assisted | `hw59_overlap.py`, `hw59_walls.py` | `python3 hw59_overlap.py`; `python3 hw59_walls.py` | `logs/hw59_overlap.txt`: distances `[2, 3]` for AA, AB, BA, BB and `bad 0`; `logs/hw59_walls.txt`: `overlap failures 0` for the lengths 5 to 9 | 1 s | quick |
| 21 | Example 7.11: the statements about frames on all binary words of length 31 and 21 | computer-assisted | `marker_S3.c` | `gcc -O2 -o marker_S3 marker_S3.c && ./marker_S3` | `logs/marker_S3.txt` | 11 min | full |
| 22 | After Proposition 6.24: the variant of Ψ_2 of radius 6 for 219R commutes with F on all 4^15 windows and is an involution | computer-assisted | `rot_psi2.c` | `gcc -O2 -fopenmp -o rot_psi2 rot_psi2.c && ./rot_psi2` | `logs/rot_psi2.txt` | 4 min on one core | full |
| 23 | Remark 7.5: Λ_M is strongly connected for M ≤ 4 for the affine rules 15, 60, 105, 150, 195, 240, and not for 90 and 165; for rule 90, for M ≤ 8 | computer-assisted | `window_graph.py` | `python3 window_graph.py`; `python3 window_graph.py 90 8` | `logs/window_graph.txt`: one component for 240, 15, 60, 195, 150, 105 and M ≤ 5; `logs/window_graph_90.txt`: `[1, 2, 4, 8, 16, 32, 64, 128, 256]` | 5 s; 15 s | quick |
| 24 | After Theorem 1.4: the orders o_L of F_L (L ≤ 7) and the periods of the content 0^L (L ≤ 10) | computer-assisted | `walls54.py`, `cycles54.py` | `python3 walls54.py`; `python3 cycles54.py` | `logs/walls54.txt`: orders 3, 3, 6, 90, 2520, 102702600, 4354384834800; periods 3, 3, 6, 9, 18, 21, 18, 21, 30, 33; `logs/cycles54.txt` | 4 s (1 GB); 3 s | quick |
| 25 | Table 2: spreading of a changed cell, from simulation on random configurations | computer-assisted (sampled) | `walls_search.py` | `python3 walls_search.py 54 110 30` | `logs/walls_search.txt` | 2 min | full |

## Confirmations of results proved by hand (Appendix A.4)

| # | statement | status | program | command | expected output / log | time, memory | tier |
|---|---|---|---|---|---|---|---|
| 26 | Theorem 1.3 on 3000 random cyclic configurations; nontrivial on 476 | confirmation only | `sym54.py` | `python3 sym54.py` | `logs/sym54.txt`: `commutes with F: True; sigma^-1 inverts: True; homomorphism: True; nontrivial on 476/3000 samples` | 1 s | quick |
| 27 | Theorem 1.1 for each rule of Table 1 on 1500 random cyclic configurations | confirmation only | `trains_true.py` | `python3 true_classes.py && python3 trains_true.py` | `logs/trains_true.txt`: verified for 31 rules, `not verified: []` | 6 s | quick |
| 28 | Theorem 1.5: the table of τ differs from the identity on 24 of 16 384 entries, commutes with F on all 4^9 windows; τ is injective (pair graph; a second test on periodic configurations) | confirmation only | `tau.py` | `python3 tau.py` | `logs/tau.txt`: `changed entries: 24`, `commutes with F: True`, `injective (pair graph): True ; periodic<=10: True` | 30 s, 4.5 GB | full |
| 29 | Theorem 1.5: τ∘τ = id on all 4^13 windows of length 13, and τD differs from Dτ on 147 456 of them | confirmation only | `tau_exact.py` | `python3 tau_exact.py` | `logs/tau_exact.txt`: `tau o tau == id on all 13-windows: True`, `windows where they differ: 147456` | 10 s | quick |
| 30 | Proposition 5.2: D commutes with F_133 on all 4^11 windows | confirmation only | `verify133.py` | `python3 verify133.py` | `logs/verify133.txt`: `commutes with F on all 4^11 windows: True` for 133 and 161 | 30 s, 1.1 GB | quick |
| 31 | Proposition 6.16: the elements of T_R composed exactly on windows of radius 4; the eight generators commute with F exhaustively at radius 3 | confirmation only | `hw59_census.py`, `hw59_generators.py`, `hw59_classes.py` | `python3 hw59_census.py`; `python3 hw59_generators.py`; `python3 hw59_classes.py` | `logs/hw59_census.txt` (row 8); `logs/hw59_generators.txt`: `True` in the eight lines; `logs/hw59_classes.txt`: `class changes 0` | 2 min; 1 s; 1 s | full; quick; quick |
| 32 | Theorem 6.14 (59R): rings of length 2000 seeded with trains of length up to 6, all m = 2, …, 6 at once, 60 steps; radius bound 2m+2 | confirmation only | `hw59_trains.py`, `hw59_table.py` | `python3 hw59_trains.py`; `python3 hw59_table.py` | `logs/hw59_trains.txt`: `Theorem 6.14: commute failures 0 hom/inverse failures 0 nontrivial 200 of 200`, `radius 2m+2 violations 0`; `logs/hw59_table.txt` ends `flip ok` | 10 s | quick |
| 33 | Proposition 6.19 (127R, 179R, 223R): the flip of all walls composed exactly at radius 4, pair-graph test | confirmation only | `hw_flipall.py`, `hw_families.py` | `python3 hw_flipall.py 127,179,223` | `logs/hw_flipall.txt`: `commutes True tau^2=id True pairgraph injective True` for each family | 9 s | quick |
| 34 | Proposition 6.19: the cycles in all contexts; rings of length 1500 with trains at gaps 5, 6, 8 and walls of other classes mixed in | confirmation only | `hwtp_cycles.py`, `hwtp_walls.py`, `hwtp_trains.py` | `python3 hwtp_cycles.py`; `python3 hwtp_walls.py`; `python3 hwtp_trains.py` | `logs/hwtp_cycles.txt`, `logs/hwtp_walls.txt`; `logs/hwtp_trains.txt`: `commute failures 0 hom failures 0` in the nine lines | 1 s; 1 s; 10 s | quick |
| 35 | Theorem 6.20 (27R): (W1), (H1) in all 16 contexts for F and F^-1; (H2), (H4) on words of length at most 11; β and β+1 are the only maps satisfying (H3); rings with trains at gaps 3, 6, 7 | confirmation only | `hw27_walls.py`, `hw27_beta.py`, `hw27_trains.py` | `python3 hw27_walls.py`; `python3 hw27_beta.py`; `python3 hw27_trains.py` | `logs/hw27_walls.txt`: `W1,H1 True`, `bad 0`; `logs/hw27_beta.txt`: two solutions; `logs/hw27_trains.txt`: `H4 violations 0`, `fails 0` | 4 s; 1 s; 10 s | quick |
| 36 | Theorem 6.23, Proposition 6.25 (219R, 126R): the hypotheses on all contexts by two implementations; bond invariance, π-evolution, the maps φ_s and R_s | confirmation only | `rot_walls.py`, `rot_verify.py`, `rot_hidden.py` | `python3 rot_walls.py`; `python3 rot_verify.py`; `python3 rot_hidden.py` | `logs/rot_walls.txt`, `logs/rot_verify.txt`, `logs/rot_hidden.txt` | 2 s | quick |
| 37 | Propositions 6.24, 6.25: D and E commute with F on all 4^7 windows and are injective by the pair-graph test | confirmation only | `rot_DE.py` | `python3 rot_DE.py` | `logs/rot_DE.txt` | 6 s | quick |
| 38 | Theorem 6.23: Ψ_1, Ψ_2, Ψ_3 on planted rings | confirmation only | `rot_trains.py`, `rot_verify.py` | `python3 rot_trains.py` | `logs/rot_trains.txt`: `commutation failures 0` for 219 and 126; item (5) of `logs/rot_verify.txt` | 10 s | quick |
| 39 | Theorem 7.9: Φ and Φ' on planted random rings (commutation with F, group law, commutation of the two movers); negative control for 165R; Lemma 7.10 with N = 2, three words of length 19 | confirmation only | `markers.py` | `python3 markers.py` | `logs/markers.txt`: all failure counts 0 for 90 and 165; `negative control 165, no complement frames: commutation failures 296 of 300` | 10 s | quick |
| 40 | Theorem 7.9: Φ and Φ' on all cyclic configurations of periods 10 and 11 | confirmation only | `marker_cells.py` | `python3 marker_cells.py` | `logs/marker_cells.txt` | 17 min, 8 GB | full |
| 41 | Theorem 2.4(1), (2) by exhaustion for all 256 rules and n ≤ 4 | confirmation only | `reach_thmA.py` | `python3 reach_thmA.py` | `logs/reach_thmA.txt`: `rules reading a and c: 228 violations: []` | 15 s | quick |
| 42 | Lemma 2.3: the construction with random staircase data for the 240 rules reading c and n ≤ 31 | confirmation only | `reach_stair.py` | `python3 reach_stair.py` | `logs/reach_stair.txt`: `rules 240 tests 1440 Lemma 2.3 failures 0 Theorem 2.4 failures 0` | 1 s | quick |

## Statements of the appendix without a program in this release

| statement | note |
|---|---|
| Section 3: the pair-graph test "was implemented twice, in pure Python and vectorised, and the two implementations were checked to agree at R = 1; at R ≤ 2 every map found injective was also checked to be injective on all spatially periodic configurations of period at most 2R+4" | both checks are assertions inside `centralizer.py` (function `run`) and are executed by rows 1 and 7; no separate program |
| Section 4.6 (row 6): the enumeration at radius 3 for 54R | the program is included, its output is not |
| Appendix A.4, 127R/179R/223R: "rings … showed no failure of commutation with F and F^-1" at gaps 5, 6, 8 | covered by row 34 |
