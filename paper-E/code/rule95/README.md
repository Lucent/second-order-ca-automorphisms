# Rule 95: block maps commuting with F_95

Cells are c_i = (x_i, y_i), and

    F(c)_i = (y_i, x_i XOR f(y_{i-1}, y_i, y_{i+1})),   f = Wolfram rule 95.

The programs in this folder concern the automorphism h_* of radius 7, the
enumeration of the conserving block maps of the golden mean shift, the
search for block maps of sublattice type with a SAT solver, the comparison
with the census at radius 3 and 4 (`../census`), and blocking words and
window graphs.

All running times are for one core unless stated otherwise. The directory
`logs/` contains the outputs of the runs named below.

## Notation

The programs are written in Python 3 (with numpy), except `verify7.c`.
`core.py` is a module used by the others: it contains F, the block maps of
the golden mean shift, and the evaluation of Ĥ_{φ,ψ} σ^n on finite
perturbations of the background (1,1). Run on its own (`python3 core.py`)
it checks on 300 pseudo-random finite perturbations each that Ĥ_{φ,ψ} is
σ^{−1}F, σF and σ^2 for (φ,ψ) = (σ,id), (id,σ) and (σ^{−1},σ), and prints
three lines ending in `mismatches 0`. The shift is (σc)_i = c_{i+1}.

- G is the golden mean shift (binary sequences without 11). A block map of G
  with window w and anchor a is (φx)_j = f(x_{j+a} … x_{j+a+w−1}); a word is
  stored as a binary number whose most significant bit is its first letter.
  "Conserving" means: the number of 1s of every finite point is preserved.
- φ_* is the involution of G exchanging the words 0001001 and 0000101
  wherever they occur, and φ̄_* the one exchanging 1001000 and 1010000. In the
  programs they are `gm_swap('0001001', '0000101')` and
  `gm_swap('1001000', '1010000')` (variables `phi1`, `phi1m`); in printed
  output they are called `phi_*`/`phi*` and `phibar_*`/`phibar*`, and σ^k is
  written `sigma^k` or `s^k`.
- Classes. A point of Φ^+ (x_i = y_{i+1} for all i) has the diagram
  π(q) = y_{−q}, and a point of Φ^− (x_i = y_{i−1}) the diagram ξ(p) = y_p.
  The two class sequences of a point of Φ^+ are x^e_j = 1 XOR π(2j+e),
  e = 0, 1, and those of a point of Φ^− are y^e_k = 1 XOR ξ(2k+e); they are
  points of G.
- Order of the words. All words of G are written in the order of increasing
  index j (or k) of the class sequence, as in the paper. For the right-movers
  (Φ^+) this is the order of decreasing cell index, so in a configuration the
  words of φ_* appear from right to left, on every second cell.
- Ĥ_{φ,ψ} σ^n, for conserving block maps φ, ψ of G, is the block map that
  commutes with F and that, on a configuration which is a point of Φ^+ on the
  left and a point of Φ^− on the right, first applies σ^n and then φ to both
  class sequences of the right-movers and ψ to both class sequences of the
  left-movers. `core.hatH` computes it on a finite perturbation of the
  background by running F backwards until the two parts are at least 64 cells
  apart, applying φ and ψ, and running F forwards for the same time.
- h_* = σ^{−1} Ĥ_{φ_*,id} is the automorphism of radius 7. It is also
  Ĥ_{φ_*σ, σ^{−1}} σ.
- Sublattice type. Class 0 consists of the bits y_k (k even) and x_k (k odd),
  class 1 of the bits x_k (k even) and y_k (k odd). A table pair (TX, TY) of
  radius R is of sublattice type if each output bit reads the bits of one of
  the two classes:
  - type P: y'_i = TY[w], bit m+R of w being y_{i+m} for even m and x_{i+m}
    for odd m, |m| ≤ R, and x'_i = TX[w'] with the complementary 2R+1 bits
    (x_{i+m} for even m, y_{i+m} for odd m);
  - type S: the same with even and odd exchanged.

  The maps Ĥ_{φ,ψ} σ^n are of type P for even n and of type S for odd n.

## The tables `tx7.npy`, `ty7.npy`, `tx7.bin`, `ty7.bin`

The two tables of h_*, of type S and radius 7, 2^15 entries each, one byte
(0 or 1) per entry: for a configuration c and a cell i,

- y'_i = TY[w], where bit m+7 of w is y_{i+m} for odd m and x_{i+m} for even
  m, |m| ≤ 7;
- x'_i = TX[w], where bit m+7 of w is x_{i+m} for odd m and y_{i+m} for even
  m, |m| ≤ 7;

and h_*(c)_i = (x'_i, y'_i). The files `.bin` are the raw 32768 bytes, the files `.npy` the
same arrays in numpy format.

    md5sum tx7.bin ty7.bin
    5a8b4dc6a76c188f2d5ccab5b224ddc3  tx7.bin
    fda0396a76179ee07663e93e84d27166  ty7.bin

## `gen7.py` — generator of the tables of h_*

Fills the two tables from the definition of h_* = σ^{−1} Ĥ_{φ_*,id}, by
evaluating `core.hatH` on finite configurations whose two class sequences are
pieces of a de Bruijn sequence of order 15. Every sample is compared with the
entry already filled, and every entry is filled.

    python3 gen7.py             # usage: gen7.py [output directory]

writes `tx7.npy`, `ty7.npy`, `tx7.bin`, `ty7.bin` (by default into the
directory of the script; the files written are identical to the ones
supplied). Expected output:

    pass 0 unfilled: 16154 16150 inconsistent samples: 0 of 72960
    pass 1 unfilled: 0 0 inconsistent samples: 0 of 145920
    saved; ones: 16384 16384

Running time: about 5 seconds.

## `verify7.c` — commutation and involution of h_*, on cells

Uses only the two tables and the definition of F. Checks exhaustively:

- (A) h_* F = F h_* at cell 0. Each output bit of both sides is a function of
  the 17 bits of one class in the cells −8..8; all 2^17 values are enumerated
  for each class, the bits of the other class being pseudo-random;
- (B) (σ h_*)^2 = id at cell 0. Each output bit is a function of the 29 bits
  of one class in the cells −12..16; all 2^29 values for each class;
- (C) both tests on 2^26 windows in which all bits are pseudo-random.

Uses OpenMP.

    gcc -O2 -fopenmp -o verify7 verify7.c
    ./verify7 tx7.bin ty7.bin     # usage: verify7 tx7.bin ty7.bin [FLIP [NB [NC]]]

Expected output:

    (A) class 0 (y-bit): 2^17 windows, mismatches 0
    (A) class 1 (x-bit): 2^17 windows, mismatches 0
    (B) class 0 (y-bit): 2^29 windows, violations of (sigma h)^2=id: 0
    (B) class 1 (x-bit): 2^29 windows, violations of (sigma h)^2=id: 0
    (C) 2^26 fully random windows: commutation mismatches 0, involution violations 0

Running time of the full check: 19 seconds on 32 threads (about 10
core-minutes). The outputs of the full check and of the full negative
control below are `logs/verify7.log` and `logs/verify7_negative_control.log`.

Optional arguments:

- `FLIP` ≥ 0 inverts the entry TY[FLIP] before the checks (negative
  control); `FLIP` < 0 leaves the table unchanged.
- `NB` (default 29) restricts (B) to the 2^NB windows in which the bits
  NB..28 of the class are 0; `NC` (default 26) sets the number 2^NC of
  windows of (C).

The full negative control `./verify7 tx7.bin ty7.bin 777` prints

    NEGATIVE CONTROL: TY[777] flipped
    (A) class 0 (y-bit): 2^17 windows, mismatches 8
    (A) class 1 (x-bit): 2^17 windows, mismatches 4
    (B) class 0 (y-bit): 2^29 windows, violations of (sigma h)^2=id: 32768
    (B) class 1 (x-bit): 2^29 windows, violations of (sigma h)^2=id: 16384
    (C) 2^26 fully random windows: commutation mismatches 6073, involution violations 6004

A reduced run,

    OMP_NUM_THREADS=2 ./verify7 tx7.bin ty7.bin -1 22 22

prints (about 3 seconds)

    (A) class 0 (y-bit): 2^17 windows, mismatches 0
    (A) class 1 (x-bit): 2^17 windows, mismatches 0
    (B) class 0 (y-bit): 2^22 windows, violations of (sigma h)^2=id: 0
    (B) class 1 (x-bit): 2^22 windows, violations of (sigma h)^2=id: 0
    (C) 2^22 fully random windows: commutation mismatches 0, involution violations 0

and `OMP_NUM_THREADS=2 ./verify7 tx7.bin ty7.bin 777 22 22` prints the
mismatches 8 and 4 in (A), the violations 256 and 128 in (B), and
`commutation mismatches 386, involution violations 364` in (C).

## `props7.py` — dependence set, restrictions and periodic points of h_*

Needs `core.py`. From the tables it computes

- for each of the 15 sites of the window, the number of pairs of windows,
  differing at that site only, with different table entries;
- the restrictions of the table map to Φ^+ and to Φ^−, on all 55^2 = 3025
  pairs of words of length 8 of G for the two class sequences, compared with
  σ^{−1} composed with φ_* on each class (Φ^+) and with σ^{−1} (Φ^−);
- for every cyclic configuration c of period n ≤ 9 (349524 configurations),
  whether h_* F c = F h_* c and (σ h_*)^2 c = c, whether σ h_* c ≠ c, and, for
  n ≤ 7, whether h_*(c) lies outside the ⟨σ, F⟩-orbit of c (for n = 8, 9 this
  is not tested and the last number printed is 0).

      python3 props7.py           # usage: props7.py [tx7.npy ty7.npy]

Expected output (lines shortened):

    exact dependence (offset m: number of windows in which flipping the bit at m changes the output):
      y': {'y-7': 128, 'x-6': 256, 'y-5': 1344, …, 'y+5': 256, 'x+6': 128, 'y+7': 0}
      x': {'x-7': 0, 'y-6': 256, 'x-5': 512, …, 'x+5': 1024, 'y+6': 128, 'x+7': 64}
    restriction to Phi^+ = sigma^-1 o (phi_* on each class): 3025 words, 0 failures; restriction to Phi^- = sigma^-1: 0 failures
    period 1: 4 configurations, failures of hF=Fh or (sigma h)^2=id: 0; sigma h moves 0; (n<=7) h(c) outside the <sigma,F>-orbit of c: 0
    …
    period 5: 1024 configurations, failures of hF=Fh or (sigma h)^2=id: 0; sigma h moves 310; (n<=7) h(c) outside the <sigma,F>-orbit of c: 310
    period 6: 4096 configurations, failures of hF=Fh or (sigma h)^2=id: 0; sigma h moves 0; (n<=7) h(c) outside the <sigma,F>-orbit of c: 0
    period 7: 16384 configurations, failures of hF=Fh or (sigma h)^2=id: 0; sigma h moves 7168; (n<=7) h(c) outside the <sigma,F>-orbit of c: 7168
    period 8: 65536 configurations, failures of hF=Fh or (sigma h)^2=id: 0; sigma h moves 11712; (n<=7) h(c) outside the <sigma,F>-orbit of c: 0
    period 9: 262144 configurations, failures of hF=Fh or (sigma h)^2=id: 0; sigma h moves 144396; (n<=7) h(c) outside the <sigma,F>-orbit of c: 0
    witness c = ((0, 0), (0, 0), (0, 0), (0, 1), (1, 0))  h(c) = ((1, 0), (0, 0), (0, 1), (0, 0), (0, 1))  F-period 15

For the periods 1 to 4 the numbers of failures and of moved configurations
are 0. All 14 counts of the y'-table from offset −7 to +6 are positive, and
so are those of the x'-table from −6 to +7: y'_i depends exactly on y_{i+m}
for odd m in [−7, 5] and on x_{i+m} for even m in [−6, 6], and x'_i exactly
on x_{i+m} for odd m in [−5, 7] and on y_{i+m} for even m in [−6, 6]. So
h_*(c)_i depends exactly on the cells i−7, …, i+7. Running time: about 25
seconds.

## `standalone7.py` — an explicit witness for h_*

Self-contained (needs numpy). Takes the configuration c = (00 00 00 01 10)^∞
(cells written x_i y_i, spatial period 5) and checks on its F-orbit that h_*
commutes with F, that (σ h_*)^2 = id and that σ h_* moves every point, and
that h_*(c) is not of the form σ^a F^b (c). Hence h_* is not in ⟨σ, F⟩.

    python3 standalone7.py tx7.npy ty7.npy

Expected output:

    c    = ((0, 0), (0, 0), (0, 0), (0, 1), (1, 0))
    h(c) = ((1, 0), (0, 0), (0, 1), (0, 0), (0, 1))
    F-period of c: 15
    h(F(d)) == F(h(d))      on the F-orbit: True
    (sigma h)^2 (d) == d    on the F-orbit: True
    sigma h (d) != d        on the F-orbit: True
    h(c) = sigma^a F^b (c) for some a, b  : False

Running time: under 1 second.

## `enum_golden.py` — conserving block maps of G

Enumerates, for w = 2, …, WMAX, all conserving block maps
(φx)_j = f(x_j … x_{j+w−1}) of G with window w, and decides for each whether
it is injective. Conservation is decided exactly: φ is conserving if and
only if there is an integer function λ on the words of G of length w−1 with
f(e) − e_0 = λ(e_0 … e_{w−2}) − λ(e_1 … e_{w−1}) for every word e of G of
length w; the program searches over λ depth first. Injectivity is decided
exactly on the graph of pairs of words with equal images, pruned to the
vertices that lie on bi-infinite paths.

    python3 enum_golden.py 8     # usage: enum_golden.py [WMAX], default 8

Expected output:

    w=2: G-words 3, conserving maps 2 (search nodes 3), injective 2: ['sigma^1', 'sigma^0']
    w=3: G-words 5, conserving maps 3 (search nodes 6), injective 3: ['sigma^2', 'sigma^1', 'sigma^0']
    w=4: G-words 8, conserving maps 6 (search nodes 16), injective 4: ['sigma^3', 'sigma^2', 'sigma^1', 'sigma^0']
    w=5: G-words 13, conserving maps 15 (search nodes 52), injective 5: ['sigma^4', 'sigma^3', 'sigma^2', 'sigma^1', 'sigma^0']
    w=6: G-words 21, conserving maps 54 (search nodes 245), injective 6: ['sigma^5', 'sigma^4', 'sigma^3', 'sigma^2', 'sigma^1', 'sigma^0']
    w=7: G-words 34, conserving maps 476 (search nodes 3065), injective 7: ['sigma^6', 'sigma^5', 'sigma^4', 'sigma^3', 'sigma^2', 'sigma^1', 'sigma^0']
    w=8: G-words 55, conserving maps 13516 (search nodes 135319), injective 10: ['sigma^7', 'sigma^6', 'sigma^5', None, 'sigma^4', 'sigma^3', None, 'sigma^2', 'sigma^1', 'sigma^0']

followed by two lines that list, for the two injective maps of window 8 that
are not shifts, the 6 and the 28 words on which they differ from σ^4. These
two maps are σ^4 φ_* and σ^3 φ̄_* (identified by `enum_golden_linear.py`).
Running time: about 1 second for WMAX = 7, about 25 seconds for WMAX = 8.
The output for WMAX = 8 is `logs/enum_golden_8.log`.

## `enum_golden_linear.py` — the same count by linear algebra

For w ≤ WMAX (at most 7) the conserving maps are found in a second way:
conservation is written as the linear equations Σ_e n_e(z) f(e) = |z| over
all cyclic words z of G of length at most the number of words of G of length
w−1, the system is solved over the rationals, and all 0/1 values of the free
coordinates are tried. The set of maps is compared with that of
`enum_golden.py`. Then the two injective maps of window 8 that are not
shifts are identified, and for a pseudo-random sample of 150 conserving maps
of window 8 and the 10 injective ones the result of the injectivity test is
compared with a search for two cyclic points of period ≤ 13 with the same
image.

    python3 enum_golden_linear.py 7

Expected output:

    w=2: 3 cycle equations, rank 2, free 1; brute count 2; DFS count 2; equal sets: True
    w=3: 5 cycle equations, rank 3, free 2; brute count 3; DFS count 3; equal sets: True
    w=4: 11 cycle equations, rank 4, free 4; brute count 6; DFS count 6; equal sets: True
    w=5: 29 cycle equations, rank 6, free 7; brute count 15; DFS count 15; equal sets: True
    w=6: 137 cycle equations, rank 9, free 12; brute count 54; DFS count 54; equal sets: True
    w=7: 2292 cycle equations, rank 14, free 20; brute count 476; DFS count 476; equal sets: True
    phi_* window 8 anchor -4 ; phibar_* window 8 anchor -3
    injective window-8 maps: ['other/shift', 'other/shift', 'other/shift', 'sigma^4 phi_*', 'other/shift', 'other/shift', 'sigma^3 phibar_*', 'other/shift', 'other/shift', 'other/shift']
    sample of 160 window-8 maps: injective-by-graph but colliding on periodic points: 0; non-injective with an explicit collision of period <=13: 149 of 149 non-injective

Running time: about 30 seconds.

## `inj_witness.py` — collisions for the non-injective maps

For every conserving map of window 7 and 8 that is not injective by the test
of `enum_golden.py`, it finds two distinct cyclic points of G of period at
most PMAX with the same image.

    python3 inj_witness.py 10    # usage: inj_witness.py [PMAX], default 14

Expected output:

    window 7: 476 maps, injective by the graph test 7; non-injective 469, with an explicit collision of periodic points 469 (largest period needed 9); injective maps with a collision: 0
    window 8: 13516 maps, injective by the graph test 10; non-injective 13506, with an explicit collision of periodic points 13506 (largest period needed 10); injective maps with a collision: 0

Running time: about 25 seconds.

## `sat_sublattice.py` — block maps of sublattice type, with a SAT solver

Needs the Python package `python-sat` (PySAT; `pip install python-sat`); the
default solver is CaDiCaL 1.5.3 (`cadical153`), and another solver of the
package can be named as the third argument. Also needs `enum_golden.py`.

The unknowns are the 2 · 2^(2R+1) entries of a table pair (TX, TY) of radius
R, of type P or S. The condition hF = Fh at cell 0 is written as clauses for
all values of the 2R+3 bits (cells −R−1..R+1) on which each output bit of
both sides depends. The program does not use the enumeration on G or the
maps Ĥ_{φ,ψ}; that every automorphism of X_{F_95} is of sublattice type is a
theorem of the paper.

    python3 sat_sublattice.py R all      # usage: sat_sublattice.py R all|aut [SOLVER]
    python3 sat_sublattice.py R aut

- Mode `all` enumerates all solutions of both types and counts those that
  are σ^a F^b (`named`), constant, or other.
- Mode `aut` enumerates, by blocking clauses on the variables concerned, the
  restrictions to Φ^+ and to Φ^− that occur among the solutions, keeps the
  injective ones, and then enumerates all solutions whose restrictions to
  Φ^+ and Φ^− are both injective. Every automorphism of radius at most R is
  among them.

The result is written to `sat_R<R>_<mode>.pkl` in the current directory.

Mode `all`, numbers of solutions (the constant map (1,1) is of both types):

| R | type P | type S | P + S − 1 |
|---|---|---|---|
| 1 | 2 | 5 | 6 |
| 2 | 10 | 5 | 14 |
| 3 | 10 | 37 | 46 |
| 4 | 226 | 37 | 262 |
| 5 | 226 | 2917 | 3142 |

For R ≤ 4 the last column is the number of all block maps of radius R
commuting with F_95 (see `census_cmp.py`). Expected output of
`python3 sat_sublattice.py 5 all` (times in brackets vary):

    R=5 type P: 4096 variables, 4608 equalities, 8192 and-constraints, 74406 clauses [0.1s]
      all solutions: 226; {'named': 25, 'const': 1, 'other': 200, 'other_injP': 50, 'other_injM': 50, 'other_inj_both': 0}; named expected 25; distinct (Phi+,Phi-) restrictions 226, distinct Phi+ 16, distinct Phi- 16 [1.2s]
    R=5 type S: 4096 variables, 2048 equalities, 8192 and-constraints, 69286 clauses [0.1s]
      all solutions: 2917; {'named': 36, 'const': 1, 'other': 2880, 'other_injP': 288, 'other_injM': 288, 'other_inj_both': 0}; named expected 36; distinct (Phi+,Phi-) restrictions 2917, distinct Phi+ 55, distinct Phi- 55 [13.0s]

Mode `aut`, numbers of solutions with injective restrictions to Φ^+ and Φ^−:

| R | type P | type S | total | of which σ^a F^b | others |
|---|---|---|---|---|---|
| 1 | 1 | 4 | 5 | 5 | 0 |
| 2 | 9 | 4 | 13 | 13 | 0 |
| 3 | 9 | 16 | 25 | 25 | 0 |
| 4 | 25 | 16 | 41 | 41 | 0 |
| 5 | 25 | 36 | 61 | 61 | 0 |
| 6 | 49 | 36 | 85 | 85 | 0 |
| 7 | 49 | 100 | 149 | 113 | 36 |

Expected output of `python3 sat_sublattice.py 7 aut`, without the progress
lines `... 2000 restrictions [121s]`:

    R=7 type P: 65536 variables, 73728 equalities, 131072 and-constraints, 1195142 clauses [1.5s]
      restrictions to Phi^+ among the solutions: 477 (projection on 110 variables); injective: 7 [40.3s]
      restrictions to Phi^- among the solutions: 477 (projection on 110 variables); injective: 7 [72.5s]
      solutions with injective restrictions to Phi^+ and Phi^-: 49, of which sigma^a F^b: 49 (expected named of this type: 49), others: 0 [76.3s]
    R=7 type S: 65536 variables, 32768 equalities, 131072 and-constraints, 1113222 clauses [1.4s]
      restrictions to Phi^+ among the solutions: 13517 (projection on 110 variables); injective: 10 [645.0s]
      restrictions to Phi^- among the solutions: 13517 (projection on 110 variables); injective: 10 [1334.7s]
      solutions with injective restrictions to Phi^+ and Phi^-: 100, of which sigma^a F^b: 64 (expected named of this type: 64), others: 36 [1343.4s]

The numbers 477 = 476 + 1 and 13517 = 13516 + 1 are the numbers of
conserving maps of window 7 and 8 found by `enum_golden.py`, plus the zero
map of the constant solution; 7 and 10 are the numbers of injective ones.
For R = 6 the corresponding lines read 477, 477, 49 (type P) and 55, 55, 36
(type S), with `others: 0`.

Running times, one core: under 2 seconds for `aut` with R ≤ 5 and for `all`
with R ≤ 4; about 15 seconds for `5 all`; about 30 seconds for `6 aut`;
between 25 and 55 minutes for `7 aut` (1.2 million clauses). Mode `all` was
not run for R ≥ 6.

Logs (in `logs/`):

| run | CaDiCaL 1.5.3 | MiniSat 2.2 (third argument `minisat22`) |
|---|---|---|
| `all`, R = 1, …, 5 | `sat_sublattice_all_R1-5.log` | `sat_sublattice_all_R1-5_minisat22.log` |
| `aut`, R = 1, …, 6 | `sat_sublattice_aut_R1-6.log` | `sat_sublattice_aut_R1-6_minisat22.log` |
| `aut`, R = 7 (the run behind "exactly 36") | `sat_sublattice_aut_R7.log` | `sat_sublattice_aut_R7_minisat22.log` |

The two solvers give the same numbers in every line, and the result files
`sat_R<R>_<mode>.pkl` of the two solvers contain the same sets of
solutions for every run of the table. The statement that an enumeration is
complete rests on a final UNSAT answer of the solver; no proof log was
produced for these answers, so they are not certified independently of the
solvers.

### Rule 175 (and other rules that do not read the middle argument)

    python3 sat_sublattice.py R all cadical153 175     # usage: sat_sublattice.py R all [SOLVER [RULE]]

With a fourth argument RULE, an elementary rule that does not read its middle
argument, the relation between four entries is TY[r0] = TX[r1] xor
f(TY[r2], ., TY[r3]). Only the mode `all` is available; every solution is
classified as σ^a F^b, constant or other and rechecked against hF = Fh on 20
pseudo-random periodic configurations, and the result is written to
`sat_rule<RULE>_R<R>_all.pkl`. Without the fourth argument the program is
unchanged. The file `logs/sat_sublattice_rule175.log` is the output for
R = 1, …, 7, and `logs/sat_sublattice_rule175_minisat22.log` the output of
the same runs with MiniSat 2.2 (`python3 sat_sublattice.py R all minisat22 175`). For every R it ends with a line

    rule 175 R=7: both types together {'named': 113, 'const': 0, 'other': 0}; 2R^2+2R+1 = 113

with `named` = 5, 13, 25, 41, 61, 85, 113: the solutions are the maps
σ^a F^b, there is no constant solution and no other one. The instance for
R = 7 has 65536 variables and 1196016 (type P) and 1114096 (type S) clauses;
running time about 15 seconds per type. As a control,
`python3 sat_sublattice.py 3 all cadical153 250` prints
`{'named': 25, 'const': 2, 'other': 20}`, the numbers of rule 95 (the constant
map is counted once in each type); the output is
`logs/sat_sublattice_rule250_R3.log`.

`sublattice_any.py` is a second encoder for the same search (needs
`python-sat` and numpy only):

    python3 sublattice_any.py 175 7      # usage: sublattice_any.py RULE Rmax [cap [SOLVER]]

It prints one line per radius and type with the numbers of clauses, of
solutions, of maps σ^a F^b of the type, of constants and of others; the
output for rule 175 is `logs/sublattice_any_rule175.log`, and with MiniSat 2.2
(`python3 sublattice_any.py 175 7 20000 minisat22`)
`logs/sublattice_any_rule175_minisat22.log`. The numbers of clauses
and of solutions agree with those of `sat_sublattice.py` for every R ≤ 7
(solutions 1+4, 9+4, 9+16, 25+16, 25+36, 49+36, 49+64; `constants=0
others=0`).

## `sat_post.py` — class maps in the SAT solutions

Reads `sat_R<R>_<mode>.pkl` from the current directory. For mode `all` it
checks that every class map occurring in the restriction of a solution to
Φ^+ or Φ^− is a conserving block map of G, or the zero map (which occurs for
the constant solution), and prints the sites on which these maps depend. For
mode `aut` it does the same for the restrictions enumerated, describes the
solutions that are not σ^a F^b, and writes their tables to
`sat_others_R<R>_P.npz` and `sat_others_R<R>_S.npz` (array `tables`, one row
per solution: the entries of TX followed by those of TY).

    python3 sat_post.py R all|aut

Expected output of `python3 sat_post.py 4 all`:

    R=4 type P: 226 solutions; distinct class maps occurring in restrictions: 16, conserving 15, zero map 1; minimal windows: ['(0, 0)', '(0, 3)', '(0, 4)', '(1, 1)', '(1, 4)', '(2, 2)', '(3, 3)', '(4, 4)', 'None']
    R=4 type S: 37 solutions; distinct class maps occurring in restrictions: 10, conserving 9, zero map 1; minimal windows: ['(0, 0)', '(0, 3)', '(1, 1)', '(1, 4)', '(2, 2)', '(3, 3)', '(4, 4)', 'None']

and of `python3 sat_post.py 5 aut`:

    R=5 type P Phi^+: 16 restrictions; distinct class maps 25, conserving 24, zero 1
    R=5 type P Phi^-: 16 restrictions; distinct class maps 25, conserving 24, zero 1
       found 25, named 25, others 0
    R=5 type S Phi^+: 55 restrictions; distinct class maps 55, conserving 54, zero 1
    R=5 type S Phi^-: 55 restrictions; distinct class maps 55, conserving 54, zero 1
       found 36, named 36, others 0

In every case the number of conserving maps is the number of distinct class
maps minus 1. For `7 aut` the lines `found …` read `found 49, named 49,
others 0` and `found 100, named 64, others 36`, followed by one line for each
of the 36 solutions (`logs/sat_post_R7_aut.log`). The file `sat_others_R7_S.npz` supplied here contains
these 36 table pairs (36 × 65536 bytes); it is used by `verify36.py`.
Running time: seconds.

## `verify_sub.py` — the SAT solutions are maps Ĥ_{φ,ψ} σ^n

Needs `core.py`, `enum_golden.py`, `gen7.py`, `census_cmp.py`. For every
non-constant solution in an output file of `sat_sublattice.py` it computes
the four class maps of the restrictions to Φ^+ and Φ^−, checks that they
have the form of the restrictions of Ĥ_{φ,ψ} σ^n (n = 0 for type P, n = 1
for type S) with conserving φ and ψ, generates the table of Ĥ_{φ,ψ} σ^n from
the definition as in `gen7.py`, and compares it with the solution in all
entries.

    python3 verify_sub.py 5 sat_R5_all.pkl     # usage: verify_sub.py R FILE

Expected output:

    R=5 type P: {'n': 225, 'const': 1, 'form': 225, 'conserving': 225, 'inj': 25, 'table_equal': 225}
    R=5 type S: {'n': 2916, 'const': 1, 'form': 2916, 'conserving': 2916, 'inj': 36, 'table_equal': 2916}

So the 3141 non-constant solutions of sublattice type and radius 5 are maps
Ĥ_{φ,ψ} σ^n; 61 of them have φ and ψ injective, and in the other 3080 φ or ψ
is not injective (`logs/verify_sub_R5.log`). Running time: about 15 minutes
(0.3 seconds per table).
`python3 verify_sub.py 3 sat_R3_all.pkl` prints (about 5 seconds)

    R=3 type P: {'n': 9, 'const': 1, 'form': 9, 'conserving': 9, 'inj': 9, 'table_equal': 9}
    R=3 type S: {'n': 36, 'const': 1, 'form': 36, 'conserving': 36, 'inj': 16, 'table_equal': 36}

With a file of mode `aut` the program treats the solutions that are not
σ^a F^b and prints φ and ψ for each.

## `gen36.py`, `verify36.py`, `tables36.npz` — the 36 automorphisms of radius 7

`gen36.py` takes the maps Ĥ_{φ,ψ} σ with φ ∈ {σ^{−3}, …, σ^4, φ_* σ, φ̄_*} and
ψ ∈ {σ^{−4}, …, σ^3, φ_*, φ̄_* σ^{−1}}, φ and ψ not both shifts (36 pairs),
and fills for each a table pair of type S and radius 7 from the definition,
as in `gen7.py`; a pair (φ, ψ) is kept if no sample contradicts an entry
already filled. The tables are written to `tables36.npz`, under the keys
`"<φ> | <ψ>"` (for instance `"phi* s | s^-1"`, which is h_*), each value
being the 2^15 entries of TX followed by those of TY.

    python3 gen36.py            # usage: gen36.py [outfile]

Expected output: 36 lines of the form `phi* s | s^-1 | inconsistent samples: 0`
and a last line `consistent with radius 7: […]` listing all 36 pairs.
Running time: about 2 minutes. The file `tables36.npz` supplied here is the
output of this program.

`verify36.py` (needs `core.py`, `sat_sublattice.py`, `verify_sub.py`,
`census_cmp.py`, `gen36.py`; the package `python-sat` is not needed) checks
each of the 36 table pairs exactly: (A) hF = Fh on all 2 · 2^17 windows of
one class, for both fillings of the other class; (B) the class maps of the
restrictions to Φ^+ and Φ^−, on all words, are those of Ĥ_{φ,ψ} σ; (C) the
set of offsets on which the tables depend. A block map commuting with F_95
is determined by its restrictions to Φ^+ and Φ^−, so (A) and (B) show that
the table map is Ĥ_{φ,ψ} σ. (D) If `sat_others_R7_S.npz` is present, the set
of the 36 table pairs is compared with the set found by the SAT search.

    python3 verify36.py         # usage: verify36.py [tables36.npz [sat_others_R7_S.npz]]

Expected output: 36 lines of the form

    phi* s   | s^-1         | commutation mismatches 0, restrictions as predicted: True, neighbourhood [-7,7]

(every line has 0 mismatches, `True` and `[-7,7]`), followed by

    tables verified: 36 all ok: True ; pairwise distinct: True
    SAT search: 36 automorphisms of sublattice type S and radius <= 7 outside <sigma,F>; equal to the set of the 36 tables above: True

Running time: about 3 seconds.

## `census_cmp.py`, `census_R3.json` — the censuses at radius 3 and 4

`census_R3.json` contains all block maps of radius R ≤ 3 that commute with
F_95 (6, 14 and 46 maps for R = 1, 2, 3): a list with one entry per radius,
each with the fields `R` and `solutions`; a solution has a `label` (`const
(1,1)`, `sigma^a F^b` or `OTHER`) and a `table`. A table is a string of
4^(2R+1) digits 0..3: the digit at position e (from 0) is the image cell
2x' + y' of the window whose cells c_{−R}, …, c_R, each written 2x + y, are
the digits of e in base 4, c_{−R} being the most significant.

The table of the census at radius 4 (262 maps) is not included because of
its size (68681990 bytes; md5 of the sorted file
`4c3163f31761c2860b476e92f83a1d53`, the order of the lines depending on the
solver).
`python3 ../census/enum_r4.py 95` writes these tables again, as the file
`raw_solutions_95_R4.txt`.
`census_cmp.py 4 FILE` expects a text file with one table per line, each a
string of 4^9 = 262144 digits 0..3 in the format just described; without
`FILE` it reads `census_R4.txt` in the directory of the script.

For every table of the census the program (needs `core.py`,
`enum_golden.py`, `gen7.py`)

1. rechecks hF = Fh on all 4^(2R+3) windows;
2. checks on all entries that the table is of sublattice type (P or S), or
   constant;
3. computes the four class maps of the restrictions to Φ^+ and Φ^−, on all
   words, and checks that they have the form of the restrictions of
   Ĥ_{φ,ψ} σ^n, n = 0 (type P) or n = 1 (type S), that φ and ψ are
   conserving and map G into G, and whether φ and ψ are injective;
4. generates the table of Ĥ_{φ,ψ} σ^n from the definition, as in `gen7.py`,
   and compares it with the census table in all entries.

       python3 census_cmp.py 3      # usage: census_cmp.py R [FILE]

Expected output: 20 lines, one for each solution in which φ or ψ is not
injective, of the form

       OTHER S phi: w=4 a=-1 tab=01000011 inj False | psi: w=4 a=-2 tab=00000111 inj True

(`tab` lists f on the words of G of length w in increasing order), and

    R=3: {'total': 46, 'commute': 46, 'const': 1, 'sublattice': 45, 'form': 45, 'conserving': 45, 'both_injective': 25, 'noninjective': 20, 'table_equal': 45, 'P': 9, 'S': 36} distinct (phi,psi,n): 45

Running time: about 6 seconds. `python3 census_cmp.py 1` and
`python3 census_cmp.py 2` print

    R=1: {'total': 6, 'commute': 6, 'const': 1, 'sublattice': 5, 'form': 5, 'conserving': 5, 'both_injective': 5, 'noninjective': 0, 'table_equal': 5, 'P': 1, 'S': 4} distinct (phi,psi,n): 5
    R=2: {'total': 14, 'commute': 14, 'const': 1, 'sublattice': 13, 'form': 13, 'conserving': 13, 'both_injective': 13, 'noninjective': 0, 'table_equal': 13, 'P': 9, 'S': 4} distinct (phi,psi,n): 13

and `python3 census_cmp.py 4 FILE` prints

    R=4: {'total': 262, 'commute': 262, 'const': 1, 'sublattice': 261, 'form': 261, 'conserving': 261, 'both_injective': 41, 'noninjective': 220, 'table_equal': 261, 'P': 225, 'S': 36} distinct (phi,psi,n): 261

(about 3 minutes, about 1 GB of memory; `logs/census_cmp_R4.log`, and
`logs/census_cmp_R3.log` for R = 3). So every non-constant map of the
two censuses is a map Ĥ_{φ,ψ} σ^n: 25 resp. 41 with φ and ψ injective (the
maps σ^a F^b) and 20 resp. 220 in which φ or ψ is not injective.

## `block_test.py` — F_95 has no blocking word

Needs `core.py`. For every word u of length L ≤ LMAX over the four cell
states and every position q in u, it looks for two configurations that
contain u at the same place and whose images under some F^n, n ≥ 0, differ
at q. First pair: c is u in the background (1,1), and c' is c with one
right-moving particle added far to the left (y_{−a} = 0, x_{−a−1} = 0), for
two positions −a of different parity. Second pair: two pseudo-random
contexts of u. For each L one line gives the number of pairs (u, q) without
a witness and the largest first time n of a difference.

    python3 block_test.py        # usage: block_test.py [LMAX], default 6

Expected output:

    length 1: 4 words x 1 positions; defect witness missing for 0 (word,position) pairs, largest first time of difference 18; random pair of contexts: missing 0, largest first time 4
    length 2: 16 words x 2 positions; defect witness missing for 0 (word,position) pairs, largest first time of difference 23; random pair of contexts: missing 0, largest first time 12
    length 3: 64 words x 3 positions; defect witness missing for 0 (word,position) pairs, largest first time of difference 30; random pair of contexts: missing 0, largest first time 15
    length 4: 256 words x 4 positions; defect witness missing for 0 (word,position) pairs, largest first time of difference 35; random pair of contexts: missing 0, largest first time 15
    length 5: 1024 words x 5 positions; defect witness missing for 0 (word,position) pairs, largest first time of difference 42; random pair of contexts: missing 0, largest first time 20
    length 6: 4096 words x 6 positions; defect witness missing for 0 (word,position) pairs, largest first time of difference 47; random pair of contexts: missing 0, largest first time 29

(5460 words, 30948 pairs of a word and a position). Running time: about 20
seconds (about 1 second for `python3 block_test.py 4`).

## `window_graph.py` — window graphs

Needs numpy. For k = +1, −1, 0 and M = 0, …, MMAX it builds the graph
W^{(k,1)}_M on the 4^(M+1) words on the cells 0..M, with an edge from each
word to each of its images under σ^k F, the cells outside the window ranging
over all values, and tests whether it is strongly connected. The same is
done for rule 90 with M ≤ 4, as a control.

    python3 window_graph.py      # usage: window_graph.py [MMAX], default 7

Expected output:

    rule 95R, W^(+1,1)_M strongly connected for M=0..7: [True, True, True, True, True, True, True, True]
    rule 95R, W^(-1,1)_M strongly connected for M=0..7: [True, True, True, True, True, True, True, True]
    rule 95R, W^(+0,1)_M strongly connected for M=0..7: [True, True, True, True, True, True, True, True]
    rule 90R, W^(+1,1)_M strongly connected for M=0..4: [True, False, False, False, False]
    rule 90R, W^(-1,1)_M strongly connected for M=0..4: [True, False, False, False, False]
    rule 90R, W^(+0,1)_M strongly connected for M=0..4: [True, True, True, True, True]

Running time: about 1 second.

## Programs and statements

| statement | programs |
|---|---|
| the numbers of conserving and of injective conserving block maps of G with window w ≤ 8 | `enum_golden.py`, `enum_golden_linear.py` (w ≤ 7), `sat_sublattice.py` with `sat_post.py` (w = 7, 8) |
| collisions for the non-injective conserving maps | `inj_witness.py` |
| the automorphism h_* of radius 7 (tables, commutation, involution, restrictions, dependence, periodic configurations, witness) | `gen7.py`, `verify7.c`, `props7.py`, `standalone7.py` |
| no automorphism of sublattice type outside ⟨σ, F⟩ of radius ≤ 6; the 36 at radius 7 | `sat_sublattice.py` (mode `aut`), `sat_post.py`, `gen36.py`, `verify36.py` |
| the commuting maps of sublattice type of radius ≤ 5 are maps Ĥ_{φ,ψ} σ^n or constant | `sat_sublattice.py` (mode `all`), `verify_sub.py` |
| the commuting maps of radius ≤ 4 are maps Ĥ_{φ,ψ} σ^n or constant | `census_cmp.py` |
| F_95 has no blocking word (words of length ≤ 6) | `block_test.py` |
| window graphs W^{(±1,1)}_M, W^{(0,1)}_M for M ≤ 7 | `window_graph.py` |
| block maps of sublattice type commuting with F_175, R ≤ 7 | `sat_sublattice.py` (fourth argument 175), `sublattice_any.py`, each with CaDiCaL and with MiniSat; logs `logs/sat_sublattice_rule175*.log`, `logs/sublattice_any_rule175*.log` |
