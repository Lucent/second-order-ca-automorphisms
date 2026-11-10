# Programs for the computations of the paper

This directory contains the programs for the computations specified in the appendix
"Computations". The sections below are: the certificate (T_P) and the examples (this section and
the next three); the holonomy behind the front; the lift certificate; the SAT enumerations
(`sat/`); the rules that read two neighbours (`halfradius/`).

## The certificate (T_P) and the examples

Setting. Radius r, rule f(y_{i-r},...,y_{i+r}) = y_{i-r} XOR g(xi_0,...,xi_{2r-1}) with
xi_k = y_{i-r+1+k}, and F(c)_i = (y_i, x_i XOR f(y_{i-r},...,y_{i+r})). The shift is
sigma(c)_i = c_{i+1}.

Truth tables. A function g is encoded as an integer whose bit number
xi_0 + 2 xi_1 + 4 xi_2 + ... is the value g(xi_0, xi_1, xi_2, ...). At radius 2 this is a
4-digit hexadecimal number (for example xi_3 = 0xff00, xi_0 xi_2 = 0xa0a0), at radius 3 a
16-digit one, and at radius 1 one hexadecimal digit with xi_0 = b, xi_1 = c for
f(a,b,c) = a XOR g(b,c).

Lattices. (n,b,s) denotes the subgroup <(n,0),(s,b)> of Z^2, where (a,b) acts as
sigma^a F^b.

Requirements: Python 3 with NumPy; the programs in `sat/` and `halfradius/kernel_sat.py` also
need PySAT (`pip install python-sat`). The scripts of this section are run from this directory
with no arguments and write to standard output unless stated otherwise.

## Files

| script | what it computes | output |
|---|---|---|
| `scan_r2.py` | For each of the 65 536 functions g at radius 2, tests (T_1), (T_2), (T_3) in turn, and counts the functions by (affine?, least P). For the functions with no (T_P), P <= 3, builds Lambda_2 and Lambda_3 from the cell rule and tests strong connectivity. States are stored layer by layer as P-bit integers; sets of states are Python sets. | `out_scan_r2.txt` |
| `scan_r2b.py` | Reads from `out_scan_r2.txt` the non-affine functions with no (T_P), P <= 3, and tests (T_4) for each (and (T_5) for any that fail (T_4)). | `out_scan_r2b.txt` |
| `xcheck.py` | Compares the two (T_P) routines `scan_r2.TPtab` and `checks3.TP` on 153 functions for P = 1, 2. Counts the non-affine g, those satisfying (N), those satisfying (N) and (D), and the self-dual non-affine g, using finite differences. | `out_xcheck.txt` |
| `tp_direct.py` | (T_P) from its definition, in a different data representation: a state is stored time by time (the snapshot of all 2r layers at each t), all functions of a batch are rows of NumPy arrays, and sets of states are Boolean vectors. Run as a script it computes the least P <= 4 with (T_P) for all 65 536 functions g at radius 2, writes the full table to `out_leastP_r2.txt`, and prints the counts by (affine?, least P), the list of functions with least P = 4, and the counts of non-affine g, of those with (N), and of those with (N) and (D), computed from the algebraic normal form. | `out_tp_direct.txt`, `out_leastP_r2.txt` |
| `compare.py` | Recomputes the least P <= 4 for every g with `scan_r2.TPtab`, writes the table to `out_leastP_r2_scan.txt`, and compares it line by line with `out_leastP_r2.txt`. | `out_compare.txt`, `out_leastP_r2_scan.txt` |
| `examples.py` | For the rules of the examples (radius 1: the elementary rules 30, 45, 75, 135, 120, 180, 210, 225, 90, 150; radius 2: the named rules and ten further rules; radius 3: six rules, all given by explicit truth tables): J(g), self-duality, algebraic normal form, least P with (T_P), the stabiliser lattices of all points of least period n with the exact numbers k(L), and the (I)_q witness with least n, then least index, then least s. Tests strong connectivity of Lambda_M, built from the layer equations. Prints, for rules 30, 45, 75, 135, the complete successor relation of the column states at P = 1 (and at P = 2 for 135), the sets Im_m, and a shortest certificate: a cycle of the successor relation followed by a chain of successors ending at a state with pi = 1. | `out_examples.txt` |
| `checks3.py` | Library for general radius, with g given as a Python function of the tuple xi: F on words and on periodic configurations, Lambda_M from the cell rule, (T_P), J(g), self-duality, and a search for (I)_q witnesses. Used by `scan_r2.py` and `xcheck.py`. Run as a script it prints these data for radius 1 and for the named radius-2 rules. | standard output |
| `ca.py` | Library for radius 1: F for an elementary rule, Lambda_M, and a count of strongly connected components (Tarjan). Used by `checks3.py`. | standard output |

## Order of execution

    python tp_direct.py  > out_tp_direct.txt     # also writes out_leastP_r2.txt
    python scan_r2.py    > out_scan_r2.txt
    python scan_r2b.py   > out_scan_r2b.txt      # reads out_scan_r2.txt
    python xcheck.py     > out_xcheck.txt
    python compare.py    > out_compare.txt       # reads out_leastP_r2.txt
    python examples.py   > out_examples.txt

Running times on one core: `tp_direct.py` 10 s, `scan_r2.py` 1 min, `scan_r2b.py` 1 min,
`xcheck.py` 3 s, `compare.py` 1.5 min, `examples.py` 1 s.

## Results at radius 2 (all 65 536 functions g)

| class | least P with (T_P) | number of g |
|---|---|---|
| non-affine | 1 | 57 482 |
| non-affine | 2 | 6 586 |
| non-affine | 3 | 1 392 |
| non-affine | 4 | 44 |
| affine | 1 | 30 |
| affine | none for P <= 4 | 2 (xi_3 = 0xff00 and 1 + xi_3 = 0x00ff) |

There are 65 504 non-affine g; 64 512 of them satisfy (N), and 64 288 satisfy (N) and (D).
The tables `out_leastP_r2.txt` and `out_leastP_r2_scan.txt` are identical.

## The holonomy behind the front (Appendix B of the paper)

These four programs confirm Lemma `pm:lem:theta` and Proposition `pm:prop:odd` of the appendix on
the two-copy window graph. They use the cell rule F(c)_i = (y_i, x_i XOR f(y_{i-1}, y_i, y_{i+1}))
only; the layer equation and the conditions (C_j) are not used. Here RULE is the number of an
elementary rule f(a,b,c) = a XOR g(b,c), and W, M are the parameters w, M of the window graph
Gamma_{w,M}. G_w is the eventual image of Gamma_{w,M} with the layer y[w] forgotten. The value
Psi of a closed walk is an element of (Z/2)^5, printed as an integer whose bits 0,...,4 are the
entries (T, l_a, l_b, N_a, N_b). In the programs `q`, `kcoef` (or `k`) and `alpha` stand for g,
g_11 and g_10.

| script | what it computes |
|---|---|
| `ref_check.py RULE W M` | All window states, the eventual image by iteration (its size is compared with 2^(3w+2) 4^M), strong connectivity of G_w, and the exact set of values Psi_w of the closed walks at a vertex, by a search in G_w x (Z/2)^5. |
| `ref_check.py RULE W M flip` | Tests that the two top flips and the swap map vertices of G_w to vertices and edges to edges, with the labels and the changes of psi_w stated in Lemma `pm:lem:theta`. |
| `ref_fast.py RULE W M` | The eventual image by breadth-first search from one of its states, the strongly connected components of G_w, and the group generated by the values Psi_w of closed walks, from the potentials of a spanning tree. For the larger cases (135 with w = 6, 7, M = 0 and (w,M) = (5,1), (6,1), (4,2)). |
| `ref_affine.py RULE W M` | As `ref_check.py`, with the group generated by Psi_w over the closed walks of all strongly connected components (for the affine rules 90, 105, 150). |
| `ref_base.py base` | Checks (Z1) and (Z2) for the period-3 history below 2 of rule 135, and that it is a closed walk of length 3 in G_2, with its value Psi_2, for several (w,M). |
| `ref_base.py lift W M` | For rule 135, samples closed walks of G_{w-1} with l = (0,0) and tests that they lift to closed walks of G_w of the same length. |
| `ref_base.py percomp RULE W M` | The number of values Psi_w of closed walks at sampled vertices (for rules whose graph is not strongly connected). |

Results. For rule 135, for 2 <= w <= 7 with M = 0 and for (w,M) = (2,1), (3,1), (4,1), (5,1),
(6,1), (2,2), (3,2), (4,2), (2,3): the eventual image has 2^(3w+2) 4^M states, G_w is strongly
connected, and all 32 values of Psi_w occur. The same holds for the rules 30, 45, 75, 120, 180,
210, 225 for w = 2, 3, 4, M = 0. For the affine rules with w = 2, 3, 4, M = 0, the group
generated by Psi_w over all closed walks has order 8, 8, 8 for rule 90 (4, 16, 64 components) and
4, 2, 4 for rules 105 and 150.

## The lift certificate (Proposition `pm:prop:cert`)

`lift_cert.py RULE NMAX TMAX D` computes the node sets R_0, R_1, ... of the subsection "The lift
certificate" for the elementary rule RULE with the parameters n_max = NMAX, T_max = TMAX, up to
depth D or to the first repetition R_{k0+rho} = R_{k0}. For each k it prints |R_k| and
test(R_k); each new set is compared with all earlier ones (as sets while they have at most
20 000 nodes, and through SHA-256 digests of the sorted node lists otherwise; different digests
prove that two sets are different). The last two lines state whether all tests hold and whether
a repetition occurs.

| command | expected last lines | time, one core |
|---|---|---|
| `python lift_cert.py 30 2 1 50` | `test(R_k) holds for all 0 <= k <= 2: True`; `repetition R_3 = R_2` (compared as sets). |R_k| = 5, 4, 8. | < 1 s |
| `python lift_cert.py 45 2 4 50` | `test(R_k) holds for all 0 <= k <= 22: True`; `repetition R_23 = R_19` (compared as sets). |R_0|, ..., |R_22| = 4, 6, 8, 24, 30, 78, 198, 313, 515, 1012, 1441, 1729, 1911, 2912, 3245, 2605, 2509, 3508, 3505, 2629, 2521, 3524, 3513. | < 1 s |
| `python lift_cert.py 135 4 12 105` | `test(R_k) holds for all 0 <= k <= 104: True`; `no repetition: the sets R_0, ..., R_105 are pairwise distinct`. |R_0| = 15, |R_104| = 200 512, |R_105| = 223 648; the largest set has 292 432 nodes. The full output is in `out_lift_cert_135.txt`. | 10 min |

## SAT enumerations (`sat/`)

Run from the directory `sat/`. Cells are 2x + y. A function g of xi_0, ..., xi_{2r-1} is given
as the string `gtab` of its 2^(2r) values, with xi_0 the most significant index bit (so
xi_0 xi_2 = `0000000000110011`, maj(xi_0, xi_2, xi_3) = `0001000101110111`,
xi_0 + xi_3 = `0101010110101010`, xi_1 = `0000111100001111` at radius 2). This is a different
convention from the hexadecimal truth tables above.

| script | what it computes |
|---|---|
| `centralizer.py --rules 30,45 --R 1,2` | For an elementary rule f and a radius R: all block maps of radius R on the four-letter alphabet that commute with F, by the encoding of the section "Enumeration of commuting block maps". Every solution is re-checked on all words, compared with the tables of sigma^a F^b, |a| + |b| <= R, and with the constant maps, and tested for injectivity by the pair graph. Asserts that every sigma^a F^b with |a| + |b| <= R is found. |
| `invertible_sat.py --rules 45 --R 1` | The same encoding with the balance condition added; every balanced solution is tested by the pair graph. |
| `cent_r.py ECA<n> R`, `cent_r.py <gtab> r R` | The same enumeration for a left-permutive rule of radius r. Prints one JSON record: `nsol` (number of commuting block maps), `named_found` (the pairs (a,b) of the solutions equal to sigma^a F^b; the comparison is exact), `extra_automorphisms` (injective solutions not of this form), `non_injective`, `constant`. Uses `common.py`. |
| `strip.py <gtab> r m` | The nonzero functions gamma of theta_0, ..., theta_m that satisfy the stripping identity (gen:eq:star), as a SAT problem on the truth table of gamma; the enumeration stops after 20 solutions. Prints their number and whether each is the constant 1. |

Expected output.

- `centralizer.py --rules 30,45 --R 1,2`: rule 30 has 7 and 15 solutions (5 and 13 maps
  sigma^a F^b and the two constants (0,0), (1,1)); rule 45 has 6 and 14 (5 and 13 maps
  sigma^a F^b and the constant (1,1)). For rule 30, `--R 3` gives 27 solutions (25 maps
  sigma^a F^b and the two constants; 45 s).
- `cent_r.py ECA45 R`, R = 1, 2: `nsol` 6, 14; `extra_automorphisms` 0; `non_injective` 1
  (`constant` 1). `cent_r.py ECA135 R`, R = 1, 2: `nsol` 5, 13, nothing else.
  `cent_r.py ECA75 1`: `nsol` 6, `extra_automorphisms` 0, `constant` 1.
- Column "extra" of Table `gen:tab:r2` (`extra_automorphisms` for r = 2 and R = 1, 2):
  xi_0 xi_2: 0, 0 (at R = 2 the nine named maps include sigma F and sigma F^-1, which have
  radius 2 because g does not depend on xi_3); maj(xi_0, xi_2, xi_3): 3, 7;
  xi_0 + xi_3: 0, 4; xi_1: 5 at R = 1.
- Column "strip" (`strip.py <gtab> 2 m`, m = 1, 2): xi_0 xi_2: 0, 0; maj: 1, 1 (the
  constant 1); xi_0 + xi_3: 15, 20 (the cap).

Running times on one core: R = 1 under 2 s; r = 2, R = 2 between 20 s and 2 min per rule;
`strip.py` under 10 s.

## Rules that read two neighbours (`halfradius/`)

| script | what it computes |
|---|---|
| `hr_checks.py [MW [MS]]` | For f = 187 and 119: (1) the permutations phi_0, phi_1 of Table `hr:tab:phi`; (2) the cell permutations commuting with F (Lemma `hr:lem:P`); (3) the sections at x_* of the elements t_1, t_2 of Lemma `hr:lem:sections`, by the product rule with free reduction, and a test of t(x_* w) = x_* t|_{x_*}(w) on all words w of length 7; (4) the order of <phi_0^2, phi_1> and its 2-transitivity (Lemma `hr:lem:level1`); (5) strong connectivity of Lambda^0_M and W_M for M <= MW (default 6), built from Definition `hr:def:graphs`; (6) the order of the level-1 image of { t|_u : t in G, t(u) = u }, u = x_*^M, for M <= MS (default 10), from Schreier generators. |
| `kernel_sat.py RULE A1,A2,...` | For a rule that does not read its left neighbour and each anticipation A: all maps k(c)_i = K(c_i, ..., c_{i+A}) that commute with F and are bijective in c_i, by SAT. Each solution is compared with the tables of F^b, |b| <= A, and the others are tested for injectivity on the configurations of period at most 6. |

Expected output.

- `python hr_checks.py`: for both rules, (2) only `id`; (3) both sections equal the stated
  ones and the test on words of length 7 is `True`; (4) order 24 for 187 and 12 for 119, both
  2-transitive; (5) `True` for every M <= 6; (6) orders 24 for 187 at M = 0, ..., 10, and
  24, 12, 12, ..., 12 for 119. Time 7 s. The output is in `halfradius/out_hr_checks.txt`.
- `python kernel_sat.py 187 1,2,3,4,5,6` and `python kernel_sat.py 119 1,2,3,4,5,6`: for each A,
  2A + 1 solutions, the powers F^b with |b| <= A, and no other solution. Time 5 s.
- `python kernel_sat.py 102 1,2,3`: 12, 48, 192 solutions, of which 6, 30, 114 are not
  injective.
