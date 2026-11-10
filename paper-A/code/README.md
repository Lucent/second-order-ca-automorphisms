# Code and output

This directory contains the programs for the computer checks of the paper, together with their outputs (`out_*.txt`). No computer is used in any proof of the paper: every program in this directory is a confirmation of a statement proved by hand. The confirmations for the class of 79 are in `c79/`. Requirements: a C compiler (`gcc`), Python 3, and `sh` with `bc` for the timings. The programs for the classes of 110, 14, 46, 24 and 189 use only the Python standard library; the programs in `c79/` use `numpy`; `figures.py` uses `numpy` and `matplotlib` (if they are not installed for `python3`, `run_all.sh` runs `figures.py` through `uv run --with numpy --with matplotlib`). Run everything with

    ./run_all.sh

and the confirmations for the class of 79 with `c79/run_all.sh`.

`run_all.sh` takes about six minutes on one core, `c79/run_all.sh` about seven minutes. Each program prints `ALL CHECKS PASSED` (or, in `c79/`, the final line listed in `c79/README.md`) as its last line, or prints `FAILURE` and exits with a nonzero status.

## Classes of 110, 14 and 46: confirmations (not used in proofs)

| Script | Checks | Command | Time | Output |
|---|---|---|---|---|
| `base_cases.c` (with `moves.h`), mode `cone` | The orbits of the moves at the times abs(m) <= k-1 on the light-cone window [r-3k+1, r], k = 1..4, by union-find: for 110 the classes of Theorem (F_k); for 42 the classes of (A) and the stronger form of (B) (classes indexed by the length of the run of 11 cells at the tip); for 46 the bond classes T_k(c) | see `run_all.sh` | 22 s | `out_base_cone.txt` |
| `base_cases.c`, mode `short` | Window [r-k, r] with 6 fillings to the left (zeros, then pseudo-random): 110 and 42 for k = 5..10, 46 for k = 4..8 | see `run_all.sh` | 84 s | `out_base_short.txt` |
| `checkA_constructive.py` | The proofs of Theorem (F_k) (110) and of (A), (B) (42) run as algorithms: for every configuration on [r-k-2, r], k <= 5 (with three random fillings further left), and for 300 random configurations for k = 6..8, an explicit word of moves to a canonical configuration, built step by step as in the proof (with the claimed class at level k-1 and the claimed tip asserted at each step), replayed with the move tables | `python3 checkA_constructive.py 5 8 300` | 211 s | `out_checkA.txt` |
| `checkB_probes.py` | Independent of the move tables: moves computed from probes and debris; orbits by breadth-first search for k <= 5; every step of the proofs for 110 and 42 (the configuration b~ lies in O_{k-1}(Fb), F^{-1} b~ in O_k(b), with the stated tip) and the transfer formula | `python3 checkB_probes.py 5 1` | 22 s | `out_checkB.txt` |
| `probes.py` | The moves computed from real probes and debris. Checks the separation of the outgoing kink, the move tables (110, 42, 46), the cone of the moves, and the bounce formula (U = 2e - N + Y_d(N-e-1, e)). | `python3 probes.py 1500` | 2 s | `out_probes.txt` |
| `checks14_46.py` | The conjugacy Psi F_14 = F_42 Psi and the images of eta, zeta, lambda; the collision formula for 1 <= N <= 30; the invariance of the bond bits of 46 under F and rho | `python3 checks14_46.py` | 1 s | `out_checks14_46.txt` |

## Classes of 24 and 189: confirmations (not used in proofs)

The rule 66 and its partner 189 are simulated on finite windows by `eng24.py` (each step of F or F^{-1} shrinks the trusted window by one cell on each side); configurations are compared only where they are determined.

| Script | Checks | Command | Time | Output |
|---|---|---|---|---|
| `checks24_rules.py` | The formula f(a,b,c) = [a=b](b xor c) for 66 and the formulas for 24, 189, 231; the allowed triples of 66 and 189; the identity f(a,b,c) xor f(b,c,d) = f_46(a xor b, b xor c, c xor d); on 200 random configurations, delta F = F_46 delta for 66 and 189, F_189 = nu F_66, F_189^2 = kappa F_66^2, F_189^4 = F_66^4, and kappa F = F kappa; the kinks K, L (images under F, kappa, delta, right tails); the gluing lemma on 300 random left parts; the bounce lemma (both rows of the glued state, the tail flip and the outgoing kink) on 500 random cases | `python3 checks24_rules.py` | 5 s | `out_checks24_rules.txt` |
| `checks24_moves.py` | On 400 random tailed configurations, the moves of 66 (conditional tail flips) are carried by delta to the moves of 46 computed from probes and debris of 46, and agree with the configuration on the light cone (-inf, r-abs(m)]; the fixed points of sigma^n F^n agreeing with (0,0) to the left of a window, restricted to its first 14 cells, are exactly those of (0,0)^Z and the kinks K_{p,0}, for 66 with n = 1..5 and for 189 with n = 4 | `python3 checks24_moves.py` | 3 s | `out_checks24_moves.txt` |
| `checks24_ctrl.py` | Directly for 66, on 60 random tailed configurations and k = 1, 2, 3: a composite of moves at the times abs(m) <= k-1 reaches a configuration whose y-row differs at r-k+1 and r-k+2 | `python3 checks24_ctrl.py` | 1 s | `out_checks24_ctrl.txt` |

### Conventions

- Bit d of a word stands for the column r-d. The cells to the right of the edge r are (0,0).
- A state is the pair of rows (t-1, t). The state at time 0 is the configuration itself, with cell (x, y) written as the integer 2x+y.
- In mode `cone`, every move is computed with four different paddings to the left of the window: zeros, ones, and two pseudo-random paddings. The results must coincide on the whole window. This is the light-cone bound, checked rather than assumed.
- `base_cases.c` also checks four properties of every move: it changes no cell outside [r-|m|, r]; it never makes the edge cell (0,0); it is undone by a move at the same time; and the components equal the stated partition.

## Class of 79

See `c79/README.md`.

## Figures

`figures.py` produces the spacetime pictures of the paper (`../tex/figs/*.pdf`, or `../figs/*.pdf` if the directory `../tex/figs` does not exist) from exact simulations, checks the formulas they illustrate, and writes the log `out_figures.txt`. It is run by `run_all.sh`.
