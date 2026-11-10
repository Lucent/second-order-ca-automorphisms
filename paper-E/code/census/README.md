# The census for rule 95

Programs and logs for the lemma "census for 95": all block maps of radius
R ≤ 4 on the alphabet of four cells that commute with F_95. The programs
need the package `python-sat` (PySAT, with CaDiCaL 1.5.3 and MiniSat 2.2)
and NumPy.

| statement of the paper | program | logs |
|---|---|---|
| R ≤ 3: 6, 14, 46 commuting block maps; 5, 13, 25 of them maps sigma^a F^b, one constant, 0, 0, 20 others | `centralizer.py` | `logs/census_95_R123_cadical153.txt`, `logs/census_95_R123_minisat22.txt` |
| R = 4: 262 commuting block maps; 41 maps sigma^a F^b, one constant, 220 others | `enum_r4.py`, `classify_r4.py` | `logs/census_95_R4_minisat22.txt`, `logs/census_95_R4_cadical153.txt` |

Every enumeration was run with both solvers. The two solvers return the
same sets of tables for every R ≤ 4 (for R = 4 the sorted files of the 262
tables are identical, md5 `4c3163f31761c2860b476e92f83a1d53`). The final
answer of each enumeration, that there is no further solution, is an
UNSAT answer of the solver; no proof log was produced for it, so it is not
certified independently of the solvers.

## `centralizer.py`

    python3 centralizer.py --rules 95 --R 1,2,3 --max 5000 --out census_R3.json
    python3 centralizer.py --rules 95 --R 1,2,3 --max 5000 --solver minisat22

Encodes the commutation with F_f for all block maps of radius R on the
alphabet of four cells, enumerates all solutions with blocking clauses (at
most `--max`; the cap is not reached), rechecks every solution on all
4^(2R+3) windows, matches it against the maps sigma^a F^b with |a|+|b| ≤ R
and the constants, and decides injectivity with the pair graph. For each
radius it prints the numbers of variables, clauses and solutions and the
counts `named`, `const`, `other`, `invertible`, followed by one line per
solution. The option `--solver` names the solver (default `cadical153`),
and `--out` writes all solutions with their tables to a JSON file. The logs
contain the two summary lines of each radius. Running time: about 10
minutes, almost all of it for the recheck and the classification at R = 3.

The further counts in the summary lines and the further fields of the JSON
file are not used in the paper: `layerwise-same`, `layer-separable`,
`cross` and `joint` say whether each of the two output bits reads only one
of the two bits of the cells (the same one, with the same or with different
local rules; the other one; or neither of these), `affine` whether the
table is affine over F_2, and `t_encode`, `t_enum`, `t_classify` are
running times in seconds.

The 46 tables of radius 3 are in `../rule95/census_R3.json`, the file
written with `--out`.

## `enum_r4.py`, `classify_r4.py`

    python3 enum_r4.py 95               # MiniSat 2.2; writes raw_solutions_95_R4.txt
    python3 enum_r4.py 95 cadical153    # the same with CaDiCaL 1.5.3

Runs the same encoding at R = 4 (524288 variables, 75486342 clauses, at
most 1000 solutions; the cap is not reached), writes the tables, rechecks
the commutation of every solution with NumPy (`classify_r4.py`) and
matches the solutions against the 41 maps sigma^a F^b and the constants;
the counts `layerwise-same`, `joint`, … are as above. It does not test
injectivity. It also writes a summary `results_R4_95.json`. Measured on one
core: with MiniSat 33 minutes and 6.6 GB of memory (enumeration 22
minutes), with CaDiCaL 47 minutes and 17.6 GB (enumeration 38 minutes).
The file of the 262 tables (69 MB) is not included;
`../rule95/census_cmp.py 4 FILE` reads it.
