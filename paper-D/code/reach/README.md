# Reach of the edge-move method

Programs for the proposition on the reach of the edge-move method and for the
remark on sampled computations. `reach1.py` needs NumPy and SciPy; `moves.py`
needs NumPy and a C compiler.

## The proposition: `reach1.py` (exact)

    python3 reach1.py classes      # the 44 rules of the 11 classes
    python3 reach1.py controls     # rules 110, 14, 93, 25, 87, 159
    python3 reach1.py RULE ...

For each rule it enumerates the spatially periodic configurations of spatial
period at most 6 that satisfy the necessary condition for a stuck right edge
at some pair of adjacent columns in every phase of the orbit, builds the
window graphs of Fix(sigma^a F^b) for 1 ≤ b ≤ 4 and 0 < |a| ≤ b, and looks
for a path from a window of such a background to a recurrent component that
contains none of its windows. Expected output: `reach1_classes.txt`, in which
every line ends with `0 with a kink`, and `reach1_controls.txt`, in which
every rule has kinks (for 110 and 14 at (a,b) = (-1,1), for 93 and 87 at
(-1,3), for 25 at (-1,1), for 159 at (-2,4)). Running time: about a minute
for a rule with many backgrounds.

Expected output of `python3 reach1.py 28 111 110`:

    rule 28: 0 backgrounds kept, 0 triples (B,a,b) tested, 0 with a kink
    rule 111: 12 backgrounds kept, 72 triples (B,a,b) tested, 0 with a kink
    rule 110: 84 backgrounds kept, 192 triples (B,a,b) tested, 96 with a kink; first: ...
    total number of triples with a kink: 96

## The proposition in its stated form: `reach1_unique.py` (exact)

    python3 reach1_unique.py classes      # the 44 rules of the 11 classes
    python3 reach1_unique.py controls     # rules 110, 14, 93, 25, 87, 159

Needs NumPy only. It keeps the same backgrounds as `reach1.py` and, for every
kept background B that lies in Fix(sigma^a F^b), checks that every window of
B of length 2b has exactly one outgoing edge in the pruned window graph, so
that B is the only point of Fix(sigma^a F^b) that agrees with B on a left
half-line. Expected output: `reach1_unique_classes.txt`, in which every line
ends with `0 with a second continuation`, and `reach1_unique_controls.txt`,
in which the numbers of triples with a second continuation are those of the
triples with a kink in `reach1_controls.txt`. Running time: under a minute.

## The remark on sampled computations: `moves.py` (sampled contexts)

    gcc -O2 -shared -fPIC -o uf.so uf.c
    sh moves_all.sh                # the ten runs of moves_log.txt
    python3 moves.py RULE [--bx BITS --by BITS] --speeds a/b --kmax K

The background is given by its x-row and y-row over one spatial period
(default: the constant (0,0)). The kinks of Fix(sigma^a F^b) with this
background on their left are glued next to the edge of configurations that
agree with the background to the right of the edge, and the glued
configuration is followed in time. A kink is kept as a probe if at k = 3 it
leaves, for some state, a configuration of the same kind changed only on the
last cells; for an empty list of probes the line `lib state` shows `[]`. A
move on the last k cells and the exterior cell is tabulated on all 4^(k+1)
states; the cells further left are drawn at random in 4 independent contexts
(fixed seed), and a state is used only if all 4 give the same result
(`ctxdis` counts the discarded ones). For each value of the exterior cell
(`ext`) the program prints the number of states, the number of orbits of the
moves (union-find), and the numbers of orbits and states from which the next
column is not read (`noread`).

`moves_log.txt` contains the ten runs:

| rule (class) | background, speed | k | result |
|---|---|---|---|
| 110 | (0,0), 1 | ≤ 7 | exterior 0: two orbits of sizes 1/4 and 3/4; otherwise one orbit |
| 46 | (0,0), 1 | ≤ 8 | 2^k orbits |
| 42 (class of 14) | (0,0), 1 | ≤ 8 | k+1 orbits |
| 87 (class of 31) | (0,0), 1/3 | ≤ 6 | as for 110 |
| 162 (class of 13) | x = 00, y = 01; 1/3 | ≤ 4 | many orbits, most of them `noread` |
| 25 | x = 01, y = 10; 1 | ≤ 4 | 3, 6, 14, 43 orbits for exterior 0 |
| 74 (class of 26) | (0,0), 1 | ≤ 4 | few moves, many orbits |
| 114 (class of 58) | (1,1), 1 | ≤ 4 | 1, 4, 9, 14 orbits for exterior 0 |
| 191 | (0,0), 2/4 | ≤ 4 | 7, 16, 48, 130 orbits for exterior 0 |
| 159 | (0,0), 2/4 | ≤ 4 | no probe is kept: no kink of speed 1/2 leaves a configuration of the same kind |

Running time: under a minute per run, except 110 at k = 7, 46 and 42 at
k = 8 and 87 at k = 6 (a few minutes each).

## Helper modules

`moves.py` uses `moves_engine.py` (gluing of kinks, tabulation of the moves),
`common.py` (rules, classes, the cell rule of F), `sim.py` (evolution of
configurations with a periodic background), `sft.py` (window graphs of
Fix(sigma^a F^b) and their paths), `bgstates.py` and `legit.py` (the states of
a background and the test that a configuration is of the required kind), and
`uf.c` (union-find). `moves_all.sh` runs the ten runs of `moves_log.txt`.
