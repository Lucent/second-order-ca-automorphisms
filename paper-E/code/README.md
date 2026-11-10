# Programs

Each directory has a README with the statement of the paper that a program
supports, the command, and the expected output.

| directory | computations |
|---|---|
| `census/` | the census of the block maps of radius at most 4 commuting with F_95 |
| `rule95/` | smallest radius, commuting maps of small radius, blocking words and window graphs for 95; the search of sublattice type for 175 |
| `rule94/` | the scattering closure for 94 |
| `rule91/` | scattering of bare defects for 91; the quantities Q_1 and Q_2 |
| `rule155/` | edge bit, right transparency, front strips and closures for 155 |
| `closure/` | window closures for 91, 94 and 155 |
| `lemmaL/` | checks of the lemmas behind Lemma L |
| `kink/` | checks of the kink rigidity theorems for 94 and 95 |
| `translation95/` | the translation property for 95 |
| `sh/` | the shortening statement (Sh) for 91 and 94 |

Solution tables and other data files of more than 5 MB are not included; the
programs that produce them are. The directories `logs/` contain the outputs
of the longer runs.

Tested with Python 3.14, NumPy 2.5, python-sat (PySAT) 1.9.dev21 with the solvers
CaDiCaL 1.5.3 (`cadical153`) and MiniSat 2.2 (`minisat22`), and gcc 15.2.
