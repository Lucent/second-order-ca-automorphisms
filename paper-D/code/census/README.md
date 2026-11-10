# The censuses of commuting block maps

Programs and logs for the censuses of the paper (the theorem on rules 110 and
30, the census theorem for the open classes, the census at radius 4, and the
control rows for 54R, 90R and 150R of the census table), and for the assembly
of the classification tables. The SAT programs need the package `python-sat`
(PySAT 1.9, with CaDiCaL 1.5.3, MiniSat 2.2 and MiniCard) and NumPy. The
programs for 90R and 150R do not call a solver, but they import
`centralizer.py`, which imports `python-sat`; `classes88.py` needs no
package, and `classify_r4.py` and `identify_r4.py` need only NumPy. The raw solution tables are not included (for radius 4 one table
has 4^9 entries and the files have 1 to 70 MB per rule); `enum_r4.py` writes
them again.

| statement of the paper | program | log |
|---|---|---|
| the 88 classes; affine and outer-permutive rules; fixed words of length ≤ 3; the 23 cell permutations | `classes88.py` | `logs/classes88.txt` |
| rules 110 and 30, R ≤ 3 | `centralizer.py` | `logs/census_110_30_R123.txt` |
| rules 110 and 30, R = 4 | `enum_r4.py`, `classify_r4.py` | `logs/census_110_R4.txt`, `logs/census_30_R4.txt` |
| census of the classes without a fixed word of length ≤ 3, R ≤ 3 | `centralizer.py` | `logs/census_R123_summary.txt` |
| census at radius 4 | `enum_r4.py`, `classify_r4.py`, `identify_r4.py` | `logs/census_R4_summary.txt`, `logs/identify_R4.txt` |
| cross-check of the censuses with other solvers; check of the final UNSAT answers at R ≤ 3 | `crosscheck.py` | `logs/crosscheck_R123.txt`, `logs/crosscheck_R4_110_30.txt`, `logs/crosscheck_R4.txt` |
| 54R, R = 1 (13 maps) and R = 2 (more than 2·10^4 maps) | `centralizer.py` | `logs/census_54_90_R1.txt`, `logs/census_54_R2_full.txt` |
| 54R, R = 2, the 23 invertible maps | `invertible_sat.py`, `analyze_inv.py` | `logs/invertible_54_R2.txt` |
| 90R, R = 1, 2: 2^5 and 2^17 maps; 10 and 26 invertible, all affine | `linear_90.py`, `analyze_90.py` | `logs/linear_90_R12.txt`, `logs/analyze_90_R2.txt` |
| 90R, R = 3: 2^65 maps | `rank90_packed.py`, `span90_R3.py`, `verify_90_structure.py` | `logs/rank_90_R3.txt`, `logs/span_90_R3.txt`, `logs/verify_90_structure.txt` |
| 150R, R ≤ 3: 2^4, 2^8, 2^12 maps | `rank90_packed.py` | `logs/rank_150_R123.txt` |

## `classes88.py`

    python3 classes88.py

Needs no package. It prints the 88 classes with their members and, for each,
whether the rule is affine or outer-permutive, the first word of length at
most 3 fixed in every context, and the nontrivial permutations of the four
cells that commute with F_f. The output is `logs/classes88.txt`; its first
line is `88 classes, 256 rules`, 31 classes have a fixed word, and 57 have
none. Running time: a few seconds.

## `centralizer.py`

    python3 centralizer.py --rules 110,30 --R 1,2,3
    python3 centralizer.py --rules RULE --R 1,2,3 --max 5000

Encodes the commutation with F_f for all block maps of radius R on the
alphabet of four cells, enumerates all solutions with blocking clauses (at
most `--max`), rechecks every solution on all 4^(2R+3) windows, matches it
against the maps sigma^a F^b with |a|+|b| ≤ R and the constants, and decides
injectivity with the pair graph. For each rule and radius it prints the
numbers of variables, clauses and solutions and the counts `named`, `const`,
`other`, `invertible`. Expected for rules 110 and 30: 7, 15, 27 solutions at
R = 1, 2, 3, of which 5, 13, 25 named and 2 constants
(`logs/census_110_30_R123.txt`).

`logs/census_R123_summary.txt` contains the two summary lines of every run
with `--max 5000` at R = 1, 2, 3: 72 rules, which represent 55 of the 57
classes without a fixed word of length at most 3 (all except the classes of
223 and 255). A line ending in `(CAPPED)` means that the cap was reached; a
run at R = 3 without a second line (rules 27, 105, 150) finished below the
cap and was not classified.

## `enum_r4.py`, `classify_r4.py`, `identify_r4.py`

    python3 enum_r4.py RULE          # writes raw_solutions_RULE_R4.txt and results_R4_RULE_fast.json
    python3 identify_r4.py RULE      # RULE = 24, 189, 63, 187; reads raw_solutions_RULE_R4.txt

`enum_r4.py` runs the same encoding at R = 4 with MiniSat 2.2 (about 5 GB of
memory, 524288 variables, about 7.5·10^7 clauses, at most 1000 solutions),
rechecks the commutation of every solution with NumPy (`classify_r4.py`) and
matches the solutions against the 41 maps sigma^a F^b and the constants. It
does not test injectivity. `logs/census_R4_summary.txt` has the summary lines
for 41 rules: the 28 open classes and the classes of 14, 24, 46, 62, 63, 79,
94, 95, 126, 133, 187, 189 and 219. `identify_r4.py` identifies the solutions
that are not named for 24 and 189 (the maps kappa sigma^a F^b) and for 63 and
187 (maps sigma^a F^b whose neighbourhood lies in [-4, 4] although
|a|+|b| > 4); expected output `logs/identify_R4.txt`, with `unidentified: 0`
in all four lines.

## `crosscheck.py`

    python3 crosscheck.py --R 1,2,3 --jobs 12 --rup
    python3 crosscheck.py --R 4 --rules 110,30 --solvers glucose4,cadical153 --jobs 2
    python3 crosscheck.py --R 4 --solvers glucose4 --jobs 3

The completeness of an enumeration rests on the last answer of the solver,
that the clauses and the blocking clauses of the solutions found are
unsatisfiable. `crosscheck.py` repeats the complete enumeration with other
solvers of PySAT (by default Glucose 4, MiniSat 2.2 and CaDiCaL 1.5.3) on
all instances of the shipped logs that finish below the cap, rechecks every
table on all windows, and compares the solution sets of the solvers with one
another as sets of tables and with the shipped log: the numbers of clauses
and solutions, the counts named / const / other, and the set itself when the
log determines it (the log lists the tables, or reports only the maps
sigma^a F^b and constants, whose tables are rebuilt). With `--rup` the
unsatisfiable formula (encoding plus one blocking clause per solution) is
given to Glucose 4 with proof logging, and the clausal proof is checked by a
forward reverse-unit-propagation checker written in Python in the same file,
which does not use a SAT solver. Each line of the logs ends with `OK`, and
the last line counts the mismatches (0). The shipped logs cover 200 instances at R ≤ 3 (72, 67 and 61 at R = 1, 2, 3;
197 with the proof check, which is skipped for the rules 27, 105 and 150 at
R = 3), the rules 110 and 30 at R = 4 with two solvers, and 39 further rules
at R = 4 with Glucose 4 (`logs/crosscheck_R4.txt` is the output of three
invocations, whose commands are recorded in it). At R = 3 an instance needs
4 to 8 minutes and 3 GB with `--rup` (more than an hour for the rules with
thousands of solutions); at R = 4 an instance needs 10 to 30 minutes and 6
to 15 GB.

## 54R, 90R, 150R

    python3 invertible_sat.py --rules 54 --R 2 --native    # MiniCard cardinality constraints
    python3 linear_90.py 1,2
    python3 analyze_90.py 2 all
    python3 rank90_packed.py 90 3
    python3 rank90_packed.py 150 1        # also 2, 3
    python3 span90_R3.py
    python3 verify_90_structure.py

`invertible_sat.py` adds the balance condition as cardinality constraints and
submits every balanced solution to the pair-graph test; expected: 23
invertible maps, 13 of them sigma^a F^b. `linear_90.py` computes the kernel
of the commutation system of 90R over GF(2) (dimension 5 and 17 for R = 1, 2)
and classifies all solutions; `rank90_packed.py` computes the kernel
dimension with bit-packed elimination (65 for 90R at R = 3; 4, 8, 12 for 150R
at R = 1, 2, 3); `span90_R3.py` and `verify_90_structure.py` check that the
affine maps and the two families of maps through the eigen-factors span the
solution space of 90R for R ≤ 3.
