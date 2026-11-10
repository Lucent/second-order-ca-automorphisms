# Data that are not in the repository

The repository contains every program, every log and all small data files.
The following files are written by the programs and are not included, because
of their size or because they are intermediate files. They will be deposited
on Zenodo with the first release. No file in the repository exceeds 5 MB.

| file | size | written by | command (run from the folder of the program) | time, memory | checksum |
|---|---|---|---|---|---|
| `paper-D/code/census/raw_solutions_RULE_R4.txt`, one per rule of the radius-4 census (41 rules) | 1 to 70 MB each | `enum_r4.py` | `python3 enum_r4.py RULE` | 5 to 35 min and about 5 GB per rule (MiniSat 2.2) | none recorded; the order of the lines depends on the solver |
| `paper-E/code/census/raw_solutions_95_R4.txt` (the 262 block maps of radius 4 commuting with 95R) | 68 681 990 bytes | `enum_r4.py` | `python3 enum_r4.py 95` (MiniSat 2.2) or `python3 enum_r4.py 95 cadical153` | 33 min, 6.6 GB; 47 min, 17.6 GB | MD5 of the sorted file: `4c3163f31761c2860b476e92f83a1d53` (the order of the lines depends on the solver) |
| `paper-E/code/rule95/sat_R<R>_<mode>.pkl`, `sat_rule<RULE>_R<R>_all.pkl` | small for R ≤ 5 | `sat_sublattice.py` | `python3 sat_sublattice.py R all` or `R aut` (optionally `SOLVER RULE`) | seconds for R ≤ 5; 25 to 55 min for `7 aut` | none recorded |
| `paper-C/code/raw_balanced_54_R3.txt` (20 001 balanced block maps of radius 3 commuting with 54R; the enumeration is stopped at this cap) | about 330 MB | `invertible_sat2.py` | `python3 invertible_sat2.py --rule 54 --R 3 --max 20001` | 19 min, 1.2 GB | none recorded; the solutions found depend on the solver |
| `paper-C/code/results/` (the incomplete enumeration at radius 3 for 54R) | not recorded | `enum54_R3.py` | `python3 enum54_R3.py run --workers 24` | not finished after days on 24 cores | none |

`paper-D/code/census/identify_r4.py RULE` reads `raw_solutions_RULE_R4.txt`
for the rules 24, 189, 63 and 187, and `paper-E/code/rule95/census_cmp.py 4 FILE`
reads `raw_solutions_95_R4.txt`.

Small tables that are included, with the checksums given in
`paper-E/code/rule95/README.md`: `paper-E/code/rule95/tx7.bin` (MD5
`5a8b4dc6a76c188f2d5ccab5b224ddc3`) and `ty7.bin` (MD5
`fda0396a76179ee07663e93e84d27166`); `make -C paper-E check` regenerates both
with `gen7.py` and compares the sums.
