# Searches for walls, eigenfields and conserved quantities

Exact finite searches behind the paragraph "Exact searches" of the section
on the open classes, and the finite checks for
the proposition "no walls in 110R". The C programs need only a C compiler;
`cq.py` needs no package, `nowalls110.py` needs NumPy.

All four searches are run on the same 38 rules: one rule of each of the 28
open classes (13, 25, 26, 28, 29, 31, 41, 43, 47, 57, 58, 61, 77, 91, 107,
109, 111, 123, 141, 143, 151, 155, 157, 159, 167, 175, 183, 191) and the
rules 14, 46, 62, 63, 79, 94, 95, 110, 133, 187. Cells are written 2x + y.

    gcc -O2 -o w1 w1.c; gcc -O2 -o hb hb.c; gcc -O2 -o factor factor.c; gcc -O2 -o conserved conserved.c

## `w1.c` — fixed and time-periodic words, invariant families (W1)

    ./w1 RULE K P S

Builds the graph on the words of length K with an edge from w to the word
that sigma^S F^P produces on the same K cells, for every value of the 2P
cells of context, and prints its strongly connected components that no edge
leaves. A proper set of words of length K that is invariant in the sense
(W1) under sigma^S F^P exists if and only if there is such a component other
than the whole graph; a word fixed by sigma^S F^P in every context is one of
size 1. The runs are (K, P) = (11, 1), (9, 2), (7, 3), (5, 4), (3, 5) with
all |S| ≤ P (`logs/w1_results.txt`, 1330 lines). Result: for 35 of the 38
rules every line reads `A^k is the only closed set`. The exceptions, all at
S = 0, are 46 and 191 for all five (K, P), and 133 for K ≥ 5
(`logs/w1_smallest_sink.txt` gives the smallest components: for 133 a
component of size 1, the wall; for 46 and 191 components of size at least 8,
the level sets of their frozen fields).

Expected output of `./w1 13 7 3 1`:

    rule 13 k 7 p 3 s 1: 1 SCCs, 1 sink SCC(s) -> A^k is the only closed set

Running time: seconds for K ≤ 9; about a minute per line for K = 11.

## `hb.c` — hidden-bit walls

    ./hb RULE K

For every position 1 ≤ j ≤ K-2 it computes the largest set of pairs (w, w')
of words of length K that differ exactly in the cell j and whose images under
F, in every context, again form such a pair; this set is nonempty if and only
if there is a family with a hidden involution at the cell j satisfying (W1)
and (H1). The runs are K = 10 (`logs/hb_results.txt`); a family with words of
length less than 10 gives one of length 10 by extension. Result: 0 surviving
pairs for every rule and every j, except rule 133 (48 pairs at j = 4 and
j = 5, the cells next to its wall 0330).

## `factor.c` — eigenfields

    ./factor RULE r p k MODE

MODE 0: the finest partition of the words of length 2r+1 compatible with
u o F^p = sigma^k o u for a block map u of radius r into any alphabet, by
union-find over all windows of length 2r+1+2p; it prints the number of
classes, and a nonconstant u exists if and only if this is at least 2.
MODE 1: the same for binary u with u o F^p = NOT sigma^k o u. The runs are
(r, p) = (5, 1), (4, 2), (3, 3), (2, 4), (1, 5) with all |k| ≤ p, both modes
(`logs/factor_results.txt`). Result: `classes 1` and `inconsistent` for all
rules except 46, 133 and 191, which have solutions at k = 0 only (frozen
fields; for 133 the indicator of the wall).

## `conserved.c`, `cq.py` — additive conserved quantities

    ./conserved RULE r
    python3 cq.py RULE r

`conserved.c` computes the rank of the linear system
h(F(w)) - h(w_1..w_r) = g(w_1..w_{r+1}) - g(w_0..w_r), w in A^{r+2}, in the
unknowns h on A^r and g on A^{r+1}, modulo the prime 32749, and prints the
dimension of the space of densities h modulo the trivial ones (`NONTRIVIAL`).
The nullity modulo a prime is at least the nullity over Q, so `NONTRIVIAL 0`
is exact; `cq.py` solves the same system over Q and prints a basis (r ≤ 3 in
reasonable time). The runs are r = 1, ..., 5 (`logs/conserved_results.txt`).
The rules with a nontrivial density, and the dimensions for r = 2, 3, 4, 5:

| rule | 26 | 46 | 77 | 91 | 94 | 95 | 123 | 133 | 175 | 191 |
|---|---|---|---|---|---|---|---|---|---|---|
| r = 2 | 1 | 1 | 1 | 2 | 1 | 1 | 2 | 0 | 0 | 0 |
| r = 3 | 1 | 2 | 3 | 3 | 1 | 2 | 2 | 0 | 0 | 1 |
| r = 4 | 1 | 4 | 3 | 3 | 1 | 4 | 2 | 1 | 1 | 2 |
| r = 5 | 1 | 8 | 3 | 3 | 1 | 4 | 2 | 1 | 1 | 4 |

All other rules of the list have `NONTRIVIAL 0` for r ≤ 5. Expected output of
`./conserved 91 2`:

    rule 91 r 2: rank 73 nullity 7 dim_h 6 trivial 4 NONTRIVIAL 2  <--

Running time: under a second for r ≤ 4, about a minute for r = 5.

## `nowalls110.py` — rule 110R

    python3 nowalls110.py

Three exhaustive checks, described in the header of the program; expected
output `logs/nowalls110.txt`: the right context first matters at t = 1, 2, 3
for 192, 48 and 16 of the 256 values of four cells; 147456 pairs whose
difference stays in its column through t = 4 and none through t = 5; and for
P = 1, ..., 6 no word whose right end cell is fixed by F^P in every right
context. Running time: about a minute.
