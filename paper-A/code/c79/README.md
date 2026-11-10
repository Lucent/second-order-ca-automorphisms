# Class of 79 (representative rule 93)

Confirmations of statements that the section on the class of 79 proves by hand (see the appendix subsection
"Confirmations for the class of 79", which specifies each one). No computer is used in the proof for this
class: every program in this folder is a confirmation. Requirements: Python 3 with numpy. `./run_all.sh` runs
everything (about seven minutes on one core) and writes the `out_*.txt` files listed below.

Libraries:

| script | content |
|---|---|
| `engine.py`, `graphs.py` | three-valued move engine and the move graphs Gamma_L on level-L states (appendix, Lemma (window)) |
| `core.py` | exact simulation of whole configurations with random left parts, and Gamma_3 built from it |

Confirmations (not used in the proof):

| script | confirms | output |
|---|---|---|
| `check_F3_hand.py` (with the argument `core`: on Gamma_3 built by `core.py`) | Lemma (the moves at the times 0,...,3), table and part (b), and the mirror rule, for every state of Gamma_3; the proof of Prop. (F*_3) followed step by step: the arrows of Steps 1-4 and their mirror images, the 64 sectors, the case distinctions of Steps 2 and 3, the table of J and the closed walks; the pairs joined by the proof are the 295,272 pairs of the proposition | `out_check_F3_hand.txt`, `out_check_F3_hand_core.txt` |
| `c_F3.py` | Prop. (F*_3) and its sharpness by breadth-first search on Gamma_3 (1152 states, 2688 edges): S(xi)\E <= R(xi) <= S(xi) for every state | `out_c_F3.txt` |
| `check_F3_second.py` (uses `engine2.py`, `graphs2.py`, `exact2.py`) | the search of `c_F3.py` repeated on Gamma_3 built by a second, independent three-valued engine and by exact simulation; the two graphs agree | `out_check_F3_second.txt` |
| `check_H2H3.py` (uses `recipe_H2.py`, `recipe_H3.py`, `recipe_states.py`) | the proofs of Props. (H_2), (H_3) run as recipes: all level-4 states for (H_2), all level-6 states (2 completed at level 7) for (H_3), and 1500 random whole configurations | `out_check_H2H3.txt` |
| `check_F4steps.py` | the transfer and Steps 1, 2, 4 of the proof of (F_4) from (F*_3), all level-4 states; closure of (F*_3) applied to b, Fb, F^-1 b gives the five sets of (F_4) | `out_check_F4steps.txt` |
| `check_bisim.py L` | greatest bisimulation G_L of Gamma_L w.r.t. the edge values (L=3: diagonal; L=2: not); with Lemma (moves commute with H) a second proof of (H_4) | `out_check_bisim_L3.txt`, `out_check_bisim_L2.txt` |
| `check_F4_graph.py k` | (F_4) directly by directed reachability in Gamma_4 (k=4 holds, k=3 fails) | `out_check_F4_graph.txt` |
| `check_kinks.py` (uses `base.py`) | vacuum, F-orbits of period-2 configurations, kinks, and the graphs of Fix(sigma^{-+1} F^3) (Lemma (kinks)) | `out_check_kinks.txt` |
| `check_gluing.py` (uses `glue.py`, `kinks.py`) | exhaustive near gluing of K^-_j to a configuration with arbitrary left part (Lemma (moves from kinks)) | `out_check_gluing.txt` |
| `check_farglue.py` | gluing on random whole configurations, forward and time-reversed (Lemma (moves from kinks)) | `out_check_farglue.txt` |
| `check_engine.py` | engine against exact simulation of the moves on random configurations | `out_check_engine.txt` |
| `check_graphs_exact.py` (uses `exact_moves.py`) | Gamma_3, Gamma_4 built by exact simulation equal the engine's graphs | `out_check_graphs_exact.txt` |
| `check_closure64.py` | 64-state closure of the relations R_l (step of Thm. (F_k)) | `out_check_closure64.txt` |

Expected final lines:
- check_F3_hand, check_F3_hand_core: `hand proof of (F*_3) confirmed on Gamma_3: True`
- c_F3: `(F*_3) holds on Gamma_3: True`
- check_F3_second: `(F*_3) holds on both graphs: True`
- check_H2H3: `(H_2) and (H_3) recipes on random whole configurations: (configuration, time) pairs 10500`
- check_F4steps: `closure: Proposition (F*_3) for b, Fb, F^-1 b joins each of the five sets of (F_4) into one class: [(((0, 0), False), 864, True), (((0, 0), True), 288, True), (((0, 1), None), 1152, True), (((1, 0), None), 1152, True), (((1, 1), None), 1152, True)]`
- check_bisim_L3: `diagonal: True pairs 1152`
- check_bisim_L2: `offdiag pairs 1372 differ at distances from edge [1, 2]`
- check_F4_graph: `F_4 states 4608 edges 13824 claimed form holds (directed, all tails): True [(1, 1152, True), (2, 1152, True), (3, 1152, True), (4, 864, True), (5, 288, True)]`
- check_kinks: `kinks ok: W lies in the big component, which has no exit; the constant loops have no exits; each vacuum pair has exactly one kink to a constant`
- check_gluing: `near-glue cases 216 branches 24896 failures 0`
- check_farglue: `rho far-glue tests 40 bad 0`
- check_engine: `z bad 0`
- check_graphs_exact: `L 4 engine graph == exact-simulation graph, 13824 edges`
- check_closure64: `(1, 1) symmetric True classes [(64, [0, 1])]`
