# Rule 155

N(a,b,c) = f(a,b,c) XOR a XOR c is 1 exactly on 000, 110, 111; the words of
Φ_155 contain neither 11 nor 000. All programs need numpy except `edge155.py`
and the C programs. `core.py` is a module (the rule, the diagram of a
configuration on a ring) used by `thmW.py`, `strip.py`, `stripsim.py` and `sonly.py`. All running times are for one core.

| statement | program | command | expected output |
|---|---|---|---|
| the edge bit: Y(t, s−t) = Y(0, s), and the 11-front fires at every step | `edge155.py` | `python3 edge155.py` | `trials 3000 violations 0 11-front cases firing every step 799` |
| right transparency (a)–(d) | `thmW.py` | `python3 thmW.py 1 3000` | `checks 99059 violations {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'phi': 0}` |
| right transparency, all fillings of w ≤ 4 cells next to the medium | `checkW.py` | `python3 checkW.py 40 4 3` | `W-check T 40 wmax 4 violations 0 checks 40920` |
| the recursion for the front strips and their eventual periodicity | `stripsim.py` | `python3 stripsim.py` | `recursion violations 0 genuine fronts 150 late strip not on a tree cycle 0` |
| periodic continuations of the 11-front | `strip.py` | `python3 strip.py 1 40` | one line per line of the strip; see below |
| periodic continuations of the 00-front at width 18 | `s4.py` | `python3 s4.py 18` | `width 18 00-front cycles 220`, `Phi- words 86 shifted 11-fronts 106 shifted 11-fronts with the final 11 on the last two lines 28 other 0` |
| window closure with base Φ^+ ∪ Φ^-, R ≤ 5 | `../closure/dc.c` | see `../closure/README.md` | `COMPLETE` |
| the same closure, R = 6 | `mkbase.py`, `lc.c` | see below | `RESULT rule 155 R 6 D 2 mode 3 unforced 0 iterations 25` (`lc155_R6.log`) |
| the closure that forces only the value at S, with and without front windows | `sonly.py`, `seeds.py`; `mkseed.py`, `sc.c` | see below | see below |

`edge155.py`: for the seeds 0..2999, a configuration that is a valid
background to the left of a cell s and pseudo-random from s on is evolved
for 60 steps.

`thmW.py SEED TRIALS`: pseudo-random configurations that agree with a point
of Φ^+ from a cell on, T ≤ 60 steps; `checks` is the number of rows examined.
`checkW.py T WMAX NMED`: for every w ≤ WMAX and every filling of w cells (both
bits) to the left of the medium, NMED pseudo-random media and left parts, T
steps each; `checks` counts the rows (1023 runs of 40 rows for the command
above).

`strip.py BETA K` builds the tree of periodic continuations line by line
from the recursion Z_m(t+1) = Z_m(t) ⊕ Z_{m−2}(t) ⊕ Z_{m−2}(t−1) ⊕
N(Z_{m−2}(t), Z_{m−1}(t), Z_m(t)), classifying each line as slaved, free or
doubling (BETA = 1: the 11-front; BETA = 0: the 00-front). Output of
`python3 strip.py 1 40` (the line numbered k is line k of the strip):

    line 1 branches 1 {'double': 1} periods [2]
    line 2 branches 1 {'double': 1} periods [4]
    line 3 branches 1 {'slaved': 1} periods [4]
    …
    line 6 branches 1 {'slaved': 1} periods [4]
    line 7 branches 2 {'free': 1} periods [4]
    line 8 branches 2 {'slaved': 2} periods [4]
    …
    line 11 branches 2 {'slaved': 1, 'double': 1} periods [4, 8]
    …
    line 39 branches 34 {'slaved': 29, 'double': 1, 'free': 2} periods [8, 16]

The window closure at R = 6:

    python3 mkbase.py 155 6 "1,1;-1,1" base6.bin      # prints: 155 6 [(1, 1), (-1, 1)] base 170
    gcc -O3 -fopenmp -o lc lc.c
    ./lc dc 155 6 2 base6.bin

`mkbase.py RULE R "a,b;…" OUT` writes the bitmap of the windows of length
2R+1 of the bi-infinite points of Fix(σ^a F^b) for the listed pairs; `lc dc
RULE R 2 BASE` iterates the same closure as `../closure/dc.c` (five stencils
on the diamond |t|+|i| ≤ 2). The log of the run is `lc155_R6.log` (about 90
minutes with OpenMP on a multi-core machine). For R = 2 and 3 the commands
with 6 replaced by R print `unforced 0 iterations 10` and `unforced 0
iterations 18` in under 2 seconds, the numbers of passes of `dc`.

The closure that forces only the value at the point S (one step back in
time) from the other four values of a stencil: `sonly.py R none` uses the
base Φ^+ ∪ Φ^- alone; `sonly.py R fronts BETAS REV SMIN SMAX` adds the
windows, centred at the offsets SMIN ≤ s ≤ SMAX from the edge, of the
periodic front strips of the fronts with edge pairs BETAS (`01`: both
kinds), and their time reversals if REV = 1. The output ends in `unforced
(n, iterations)`, n the number of windows not forced.

    python3 sonly.py 2 none                  # seed 16 unforced (968, 7)
    python3 sonly.py 3 none                  # seed 30 unforced (16140, 10)
    python3 sonly.py 2 fronts 01 1 -2 2      # seed 61 unforced (0, 15)
    python3 sonly.py 3 fronts 01 1 -3 3      # seed 128 unforced (0, 31)

For R = 4 the same closure is computed by `sc.c` from a seed file written by
`mkseed.py R BETAS REV SMIN SMAX GENUINE OUT` (BETAS `-`: no fronts):

    gcc -O3 -o sc sc.c
    python3 mkseed.py 4 01 1 -4 5 0 s4.bin; ./sc 4 s4.bin     # seed 276, RESULT R 4 seed 276 unforced 0
    python3 mkseed.py 4 - 1 0 0 0 s4n.bin; ./sc 4 s4n.bin     # seed 54, RESULT R 4 seed 54 unforced 261168

With SMAX smaller the closure is not complete: 875 windows remain for R = 2,
SMAX = 1 (`python3 sonly.py 2 fronts 01 1 -2 1`). Running times: under
10 seconds each.
