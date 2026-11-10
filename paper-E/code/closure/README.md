# The window closure

`dc.c` computes, for an elementary rule f and a radius R, the window closure
of the paper for the base Φ^+ ∪ Φ^-: the set of windows w ∈ A^(2R+1) at
which a block map h of radius R that commutes with F_f and fixes Φ^+ ∪ Φ^-
is forced to act as the identity. The base consists of the windows of the
points of Φ^+ and Φ^-; a word of Φ_f is admitted only if it extends by 6
letters on each side, which removes the words that cannot be continued (for
155 the resulting bases are those computed by `../rule155/mkbase.py`). A window is added if, for some extension by two cells on
each side, one of the five instances of the spacetime relation through the
centre forces the difference Δ = Y_{h(c)} ⊕ Y_c to vanish there, given the
points of the diamond |t| + |i| ≤ 2 whose windows are already in the set.
Each pass uses the set of the previous passes.

    gcc -O3 -o dc dc.c
    ./dc RULE R

Output: the size of the base, one line per pass, and `COMPLETE after n
passes` if every window is forced, `STALLED` otherwise. `dc_omp.c` is the
same program with an OpenMP loop (`gcc -O3 -fopenmp -o dc_omp dc_omp.c`).

Results (`logs/dc_R1-4.log`, `logs/dc155_R5.log`):

| rule | R = 1 | R = 2 | R = 3 | R = 4 | R = 5 |
|---|---|---|---|---|---|
| 91 | base 23, 2 passes | 85, 4 | 295, 8 | 1005, 12 | |
| 94 | base 22, 2 passes | 72, 4 | 226, 7 | 700, 10 | |
| 155 | base 8, 3 passes | 16, 10 | 30, 18 | 54, 25 | 96, 35 |

All are `COMPLETE`. Running times of `dc` on one core: under 5 seconds for
R ≤ 3, about 2 minutes for 155 at R = 4; R = 5 needs about 16 times the work
per pass. For 155 at R = 6 see `../rule155` (`lc.c`).
