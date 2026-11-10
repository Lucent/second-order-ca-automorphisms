# Rule 94: scattering lemmas and the scattering-closure theorem

Cells are c_i = (x_i, y_i), and

    F(c)_i = (y_i, x_i XOR f(y_{i-1}, y_i, y_{i+1})),   f = Wolfram rule 94.

The programs in this folder check, for F = F_94, the scattering identities
used in the proof that every point of Φ^+ ∪ Φ^- with finitely many kinks lies
in the scattering closure C_L of the bare points, and they compute that
closure by brute force on finite universes.

There are two sets of programs that share no code:

- `lib.py`, `cc.py`, `verify1.py`, `verify2.py`, `verify3.py`,
  `closure_bf.py`: based on the Goursat recurrence in the light-cone
  coordinates q = t − i, p = t + i, with a direct simulation of F as an
  alternative;
- `rc.py`, `ta.py`, `tc.py`, `te.py`: based on direct simulation of F only.

All programs are in Python 3; only `cc.py` needs numpy. All running times are
for one core.

## Notation

- A point of Φ^+ has a space-time diagram of the form π(q), q = t − i; a point
  of Φ^- one of the form ξ(p), p = t + i. Here π, ξ are binary sequences
  without the word 010; the programs only use sequences that are constant
  near −∞ and near +∞ ("objects").
- A kink of a sequence is a position x with different values at x − 1 and x.
  K_a is the bare kink with diagram 0 for x < a and 1 for x ≥ a; K̄_a is its
  complement; Z_a is the isolated zero (0 at a, 1 elsewhere). The bare points
  are the two backgrounds and the K_a, K̄_a.
- For a placed pair (P, X), P ∈ Φ^+ on the left and X ∈ Φ^- on the right with
  the same background between them, S(P, X) = (Γ, Ω) is the pair after the
  collision: Γ ∈ Φ^- and Ω ∈ Φ^+. A firing is a point of the diagram at which
  the nonlinear term of rule 94 is 1 (neighbours 0, centre 1).
- [b; U] is U with a kink prepended at b (b at least 2 before the first kink
  of U); ⟨d; U⟩ is U with an isolated zero prepended at d.
- C_L is the smallest set that contains the bare points and has the
  properties (C0) a point all of whose words of length L are words of points
  of C_L of the same sign lies in C_L; (C1) invariance under the shift, rev
  and μ; (C2) if P, X ∈ C_L then both components of S(P, X) are in C_L;
  (C3) if both components of S(P, X) are in C_L then P, X ∈ C_L.
- Core length. In `rc.py`, `ta.py`, `tc.py`, `te.py` the core of a sequence
  is the shortest word between its two constant tails (a bare kink has core
  length 0). In `lib.py`, `verify2.py`, `closure_bf.py` the "core" also
  contains the first letter of the right tail, so for a sequence that is not
  constant its length is larger by 1.

## First set: Goursat recurrence

### `lib.py`, `cc.py` — modules

`lib.py`: solver for the Goursat problem on a quadrant with boundary values
π and ξ (`solve`), and the scattering map on words (`scatter`).
`cc.py`: the scattering map on kink lists computed in two ways, by a direct
simulation of F in the coordinates (t, i) (`direct`, uses numpy) and by a
Goursat solver (`goursat`). Neither is run on its own.

### `verify1.py` — scattering identities

    python3 verify1.py [lib|direct|goursat_cc] [ncases]     # defaults: lib 3000

The first argument selects the implementation of S: `lib` (Goursat solver of
`lib.py`), `direct` or `goursat_cc` (the two functions of `cc.py`). It checks,
on `ncases` pseudo-random objects with 1 to 7 kinks for each of the
randomised items:

- first-kink claims: the position of the first kink of Γ for a bare kink
  probe;
- A-pass, B-pass: the two cases of the lemma on prepending a kink,
  S(bare kink, [b; U]) = ([b; Γ'], Ω') with (Γ', Ω') = S(complementary bare
  kink, U), and the bounds on the first kink of Γ';
- exchange / mirrored / zero-zero: S(K̄_a, Z_d) = (K_{d+2}, Z_{a+1}) for
  a ≡ d mod 2 with exactly two firings, the mirror statement, and the
  crossing without firing of two isolated zeros of different parity, for all
  2 ≤ a, d ≤ 11;
- seed: the seed event S(K_a, X) = (Γ, Z_{a+1}) with its two firings, for
  the gaps 2..8, each at four placements (two positions of the object and,
  for each, the two positions of the probe of the required parity): 28 cases;
- zero prepend, zero-probe pass: the two cases of the lemma on prepending an
  isolated zero.

Expected output of `python3 verify1.py`:

    implementation: lib
    first-kink claims 3000 bad 0
    A-pass 3000 bad 0
    B-pass 3000 bad 0
    exchange / mirrored / zero-zero 300 bad 0
    seed 28 bad 0
    zero prepend y>=3 3000 bad 0
    zero-probe pass 3000 bad 0

Every line ends in `bad 0` for the other implementations and other values of
`ncases` as well. Running time: about 75 seconds for the default; about
8 seconds for `verify1.py lib 300`; under 1 second for
`verify1.py direct 200`; about 3 seconds for `verify1.py goursat_cc 300`.

### `verify2.py` — the induction, run as an algorithm

    python3 verify2.py LW NMAX [lib|direct|goursat_cc]

For every object with core length ≤ NMAX (convention of `lib.py`) it runs
the proof that the object lies in C_L for L = LW: the first kink or the
first isolated zero is removed, and the chain of scatterings of the
prepending lemmas is carried out until the gap is at least LW, where (C0)
applies. Every scattering is computed and compared with the form asserted in
the lemmas; a violation raises an AssertionError.
The count `passes` is the total number of steps of these chains (at each
step two scatterings are computed, that of the object and that of the
shortened object), and `maxpass` the length of the longest chain.

Expected output of `python3 verify2.py 10 12` (about 100 seconds):

    core length 1 objects certified so far 2 {'passes': 0, 'maxpass': 0}
    core length 2 objects certified so far 3 {'passes': 0, 'maxpass': 0}
    core length 3 objects certified so far 5 {'passes': 31, 'maxpass': 16}
    core length 4 objects certified so far 9 {'passes': 84, 'maxpass': 16}
    core length 5 objects certified so far 16 {'passes': 171, 'maxpass': 16}
    core length 6 objects certified so far 28 {'passes': 354, 'maxpass': 16}
    core length 7 objects certified so far 49 {'passes': 665, 'maxpass': 16}
    core length 8 objects certified so far 86 {'passes': 1199, 'maxpass': 16}
    core length 9 objects certified so far 151 {'passes': 2135, 'maxpass': 16}
    core length 10 objects certified so far 265 {'passes': 3764, 'maxpass': 16}
    core length 11 objects certified so far 465 {'passes': 6613, 'maxpass': 16}
    core length 12 objects certified so far 816 {'passes': 11611, 'maxpass': 16}

`python3 verify2.py 8 11 direct` ends with
`core length 11 objects certified so far 465 {'passes': 4467, 'maxpass': 12}`.
A short run: `python3 verify2.py 6 8` (about 2 seconds) ends with
`core length 8 objects certified so far 86 {'passes': 451, 'maxpass': 8}`.

### `verify3.py` — closed forms of two-body diagrams

    python3 verify3.py

Checks at every point of the quadrant (40 × 40, 2 ≤ a, d ≤ 9) the closed
form of the diagram of the pair (K̄_a, Z_d), a ≡ d, with its two firings, and
of the pair (Z_a, Z_d), a ≢ d, without firing; and that the diagram of the
pair (Z_a, K̄_u), a ≡ u, equals 1 XOR [p = u+1] for q ≥ a + 2.

Expected output:

    points checked 51200 bad 0
    part (ii): points checked 20800 bad 0

Running time: under 1 second.

### `closure_bf.py` — brute-force closure

    python3 closure_bf.py L N CAP PCAP

Computes the closure of the bare points under (C0)–(C3) inside the universe
of all sequences with core length ≤ CAP (convention of `lib.py`), using the
scatterings whose probe has core length ≤ PCAP, and checks that every
sequence with core length ≤ N is reached. Sequences that leave the universe
are dropped, so the result is a lower bound for C_L. No lemma is used.

Expected output of `python3 closure_bf.py 4 9 9 2` (under 1 second):

    round 1 known 28 of 153 ; core<= 9 : 28 / 153
    round 2 known 37 of 153 ; core<= 9 : 37 / 153
    round 3 known 51 of 153 ; core<= 9 : 51 / 153
    round 4 known 82 of 153 ; core<= 9 : 82 / 153
    round 5 known 153 of 153 ; core<= 9 : 153 / 153
    missing []

`python3 closure_bf.py 6 13 13 2` (about 4 seconds) ends with

    round 10 known 1434 of 1434 ; core<= 13 : 1434 / 1434
    missing []

## Second set: direct simulation

### `rc.py` — module

Direct simulation of F_f for an arbitrary elementary rule f on a finite
array, rows stored as integers: the scattering map `scatter` with the list of
firings, the test `in_phi` that a sequence is the diagram of a point of Φ^±,
and the enumeration `objects` of all sequences with constant tails and core
length ≤ a bound. Not run on its own.

### `ta.py` — identities and inequalities of the lemmas

    python3 ta.py [RULE [NC [NRANDOM]]]       # defaults: 94 9

Checks, by simulation: that the bare kinks and isolated zeros are points of
Φ^±; the diagrams, firings and outputs of the pairs of bare kinks (crossing
without delay for different parity, delay by 1 and two firings for equal
parity); the seed event (output and exact firing set, gaps 2..12); the three
parts of the exchange lemma, with the closed form of the diagram; and, for
every object U with core length ≤ NC (or NRANDOM pseudo-random objects with
core lengths 14..40), the identities and first-kink bounds of the lemmas on
prepending a kink (distances 2, 3, 4, 5, 8) and an isolated zero (distances
3, 4, 5, 6, 9). One line is printed per statement, with the number of cases
`n` and the number of failures.

Expected output of `python3 ta.py 94 13`: first line
`rule 94 objects 2513`, then 39 lines, all with `fail=0`; among them

    lem B output                                            n=810     fail=0
    prepend c=0: S(P_i,X_i)=(K*Gamma',Omega')               n=10810   fail=0
    prepend c=1: S(P_i,X_i)=(K*Gamma',Omega')               n=14320   fail=0
    prependZ u=d+3: S=(Z_d*Gamma',Omega')                   n=2864    fail=0
    prependZ u>=d+4: S=(K_{d+2}*Gamma',Omega')              n=11456   fail=0

Running time: about 4 seconds; under 1 second for the default
`python3 ta.py` (`rule 94 objects 265`, all `fail=0`).

The line `lem B output` counts the seed event for the gaps 2..12 at all
placements −3 ≤ a, b ≤ 3 of the required parity and three corners (810
cases). The line `cor zero: X=[b;[r;K_{r+2}]] n=25` is the identity of
sequences X = [b; [r; K̄_{r+2}]] for the gap 2 at the 25 placements (a, b);
it involves no scattering beyond the seed event. The lines `pair…` and
`cross94…` are the crossings of two bare kinks for −4 ≤ a, b ≤ 4, and the
lines `exch(i)`, `exch(ii)`, `exch(iii)` the three parts of the exchange
lemma (41, 41 and 40 cases).

`python3 ta.py 94 9 400` (about 1 second) runs the checks of the two
prepending lemmas on 400 pseudo-random objects with core lengths 14..40
(first line `rule 94 objects 400`, all lines `fail=0`).

Other rules serve as controls: `python3 ta.py 90 5` prints `fail=41` for the
delayed crossing of two bare kinks and `fail=810` for the seed event, and
`python3 ta.py 95` stops with `STEP 0 FAILS` because K_a and K̄_a are not
points of Φ^± for rule 95.

### `tc.py` — the induction, run mechanically

    python3 tc.py L MAXCORE [RULE [NRANDOM]]      # default rule 94

For every object with core length ≤ MAXCORE (or NRANDOM pseudo-random
objects with core lengths 14..30) it carries out the induction step of the
scattering-closure theorem: decomposition as [b; U] or ⟨d; U⟩; the chain of
2L scatterings of the prepending lemma, each one simulated together with the
scattering of the shortened object and compared with the asserted form; the
bound on the final gap; and the window condition of (C0), position by
position.

Expected output of `python3 tc.py 10 12` (about 3 seconds):

    rule 94 L 10 maxcore 12 objects 1434 {'bare': 4, 'prepend': 1079, 'prependZ': 350, 'zero': 1} scatterings 58164 largest core seen 45 FAILURES 0

and of `python3 tc.py 12 13` (about 7 seconds):

    rule 94 L 12 maxcore 13 objects 2515 {'bare': 4, 'prepend': 1895, 'prependZ': 615, 'zero': 1} scatterings 122242 largest core seen 53 FAILURES 0

A short run, `python3 tc.py 4 8`, prints

    rule 94 L 4 maxcore 8 objects 153 {'bare': 4, 'prepend': 112, 'prependZ': 36, 'zero': 1} scatterings 2474 largest core seen 22 FAILURES 0

Further runs (each a few seconds or less):

    python3 tc.py 16 10
    rule 94 L 16 maxcore 10 objects 467 {'bare': 4, 'prepend': 349, 'prependZ': 113, 'zero': 1} scatterings 29894 largest core seen 54 FAILURES 0

`python3 tc.py L 11` for L = 1, 2, 4, 6 prints `objects 818` and
`scatterings 3824`, `7076`, `13580`, `20084`, with `FAILURES 0`.
With pseudo-random objects:

    python3 tc.py 8 0 94 1500
    rule 94 L 8 maxcore 0 objects 1500 {'prepend': 1151, 'prependZ': 349} scatterings 48982 largest core seen 71 FAILURES 0
    python3 tc.py 20 0 94 300
    rule 94 L 20 maxcore 0 objects 300 {'prepend': 226, 'prependZ': 74} scatterings 24204 largest core seen 98 FAILURES 0

### `te.py` — brute-force closure

    python3 te.py RULE NC L1,L2,... [PROPS]       # PROPS default 0123

Computes S(P, X) by simulation for every ordered pair of sequences with core
length ≤ NC and both relative parities, checks the duality
S(μΩ, μΓ) = (μX, μP) on a pseudo-random sample, and then, for each L in the
list, iterates the closure properties named in PROPS ((C0)–(C3)) to a
fixed point, starting from the bare points. Sequences that leave the universe
are dropped, so the result is a lower bound for C_L; the computation is
meaningful for L ≤ NC + 1. No lemma is used.

Expected output of `python3 te.py 94 11 4,6,8,10,12` (about 40 seconds):

    rule 94 universe: core <= 11 size 818
    scattering table: 682120 simulations 35 s
    duality S(mu Omega, mu Gamma)=(mu X, mu P): sample 1557 failures 0
    rule 94 L 4 props 0123 : closure has 818 of 818 rounds 3 first-reached-by {'C0': 80, 'C3': 154, 'C2': 580}
    rule 94 L 6 props 0123 : closure has 818 of 818 rounds 3 first-reached-by {'C0': 20, 'C3': 181, 'C2': 613}
    rule 94 L 8 props 0123 : closure has 818 of 818 rounds 3 first-reached-by {'C0': 10, 'C3': 148, 'C2': 656}
    rule 94 L 10 props 0123 : closure has 818 of 818 rounds 3 first-reached-by {'C0': 7, 'C3': 414, 'C2': 393}
    rule 94 L 12 props 0123 : closure has 818 of 818 rounds 4 first-reached-by {'C0': 2, 'C3': 199, 'C2': 613}

(the number of seconds varies). `python3 te.py 94 12 4,6,8,10,12,13` gives
`closure has 1434 of 1434` for every L (2096356 simulations, about
2.5 minutes). A short run, `python3 te.py 94 8 4,6` (about 1 second), gives
`closure has 153 of 153` for both values of L.

`python3 te.py 94 8 1,2,3,4,5,6` (about 1 second) prints
`scattering table: 23850 simulations` and `closure has 153 of 153` for
every L.

Without (C0) the closure consists of the bare points only:
`python3 te.py 94 9 4,6,8 123` gives `closure has 4 of 267` for every L.

Rule 95 as a control. With RULE = 95 the universe consists of the sequences
on the background 1 without two zeros at distance 2, and the bare points are
the background and the sequence with a single zero. The closure is then not
the whole universe: `python3 te.py 95 7 5,6` prints

    rule 95 L 5 props 0123 : closure has 9 of 16 rounds 3 first-reached-by {'C0': 5, 'C3': 2}
    rule 95 L 6 props 0123 : closure has 3 of 16 rounds 2 first-reached-by {'C0': 1}

each followed by the list of the missing sequences and by the distances
between zeros of the same class that occur in the closure.

## Programs and statements

| statement | programs |
|---|---|
| crossing of two bare kinks | `ta.py` (lines `pair…`, `cross94…`) |
| prepending a kink (two cases, first-kink bounds) | `verify1.py` (A-pass, B-pass), `ta.py` (lines `prepend c=0`, `prepend c=1`) |
| seed event | `verify1.py` (line `seed`), `ta.py` (lines `lem B…`) |
| the sequence [b; [r; K̄_{r+2}]] of the corollary on prepending a zero | `ta.py` (line `cor zero`) |
| isolated zeros (exchange lemma) | `verify1.py` (exchange / mirrored / zero-zero), `verify3.py`, `ta.py` (lines `exch…`) |
| prepending an isolated zero | `verify1.py` (zero prepend, zero-probe pass), `ta.py` (lines `prependZ…`) |
| every point with finitely many kinks lies in C_L (the induction) | `verify2.py`, `tc.py` |
| the closure of the bare points, without the lemmas | `closure_bf.py`, `te.py` |
| time-reversal duality of S | `te.py` (line `duality`) |
