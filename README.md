# second-order-ca-automorphisms

Programs, recorded outputs and paper sources for five papers on the
automorphism groups of second-order cellular automata.

For an elementary rule f, the second-order ("Fredkin") automaton

    F(x, y)_i = ( y_i ,  x_i XOR f(y_{i-1}, y_i, y_{i+1}) )

on the alphabet {0,1}^2 is reversible. The papers study the group of the
invertible cellular automata that commute with F, which is the automorphism
group of the spacetime subshift of F. A rule f is written "f R" (for example
110R) when the second-order automaton is meant. Up to reflection and the
cellwise complement there are 88 such automata.

## The five papers

Author of all five: Michael Dayah. arXiv identifiers: to be added.

**A. Reversible rule 110 has no hidden symmetries: six families of
second-order cellular automata settled by edge moves** (`paper-A/`).
The automorphism group is ⟨σ, F⟩ ≅ Z² for the fourteen rules in the classes
of 110, 14, 46 and 79, and ⟨σ, F⟩ × Z/2 for the four rules in the classes of
24 and 189. The proofs are by hand: kinks glued to a configuration with a
constant right tail act as moves of its edge cell, and a controllability
theorem for these moves determines a commuting block map one column at a
time. Every program in `paper-A/code/` is a confirmation of a statement
proved by hand. arXiv: to be added.

**B. Fronts at the speed of light: rigidity of permutive second-order
cellular automata** (`paper-B/`).
The automorphism group is ⟨σ, F⟩ ≅ Z² for rule 30, for all 16 elementary
rules that are permutive in an outer variable and not affine, and for the 8
rules in the classes of 63 and 187; for automata of any radius that are
permutive in the leftmost variable, four checkable conditions imply the same
conclusion. The computations are finite certificates for these conditions
(among them a scan of all 65 536 rules of radius 2), SAT enumerations of
commuting block maps, and confirmations. arXiv: to be added.

**C. Every finite group from a wall: large automorphism groups of
second-order cellular automata** (`paper-C/`).
A word that the automaton maps to itself in every context, with a degree of
freedom next to it, embeds every finite group in the automorphism group. This
covers 42 of the 88 classes, through fixed words, walls with a hidden bit,
time-periodic and rotating walls, and, for the affine rules 90R and 165R,
marker permutations; the paper also treats the block maps of small radius
that commute with 54R, 90R and 150R, and the affine classes. The
computations are censuses of commuting block maps, searches for walls, and
confirmations. arXiv: to be added.

**D. The 88 second-order elementary cellular automata: a classification by
automorphism group** (`paper-D/`).
Each of the 88 classes is assigned a status from the results of the other
four papers; 28 classes remain open. An exhaustive SAT enumeration, repeated
with several solvers, shows that in each open class every automorphism of
radius at most 3 lies in ⟨σ, F⟩, and that for rules 110 and 30 the block maps
of radius at most 4 commuting with F are the maps σ^a F^b and two constants.
arXiv: to be added.

**E. The cellular automata that commute with Rule 54, and rigidity of the
second-order automaton 94R** (`paper-E/`).
The automaton 95R consists of two interleaved copies of the reversible
automaton known as Rule 54; the group of the invertible cellular automata
commuting with Rule 54, and the automorphism group of 95R, are determined,
and contain every finite group. For 94R the group is ⟨σ, F⟩ ≅ Z². Partial
results are given for 91R and 155R. The computations include SAT
enumerations up to radius 7 and exhaustive closure searches.
arXiv: to be added.

## Requirements

- Python 3.12 with the packages pinned in `environment/requirements.txt`
  (NumPy; python-sat for the SAT enumerations; SciPy and Matplotlib for one
  program each);
- a C compiler (`gcc`) and `make`;
- for `make pdf`, a LaTeX installation with `pdflatex`.

`environment/README.md` gives the exact versions, and
`environment/Dockerfile` an image that contains them.

## Quick start

    pip install -r environment/requirements.txt
    make check

`make check` runs the quick checks of every paper (`make -C paper-X check`
for a single paper) and ends with the line `all quick checks passed`. It
works in scratch copies (`paper-X/.check/`) and does not modify the committed
files. If python-sat is not installed, the steps that need a SAT solver are
skipped with a message and the others still run. To use another interpreter:
`make check PYTHON=/path/to/python`.

    make pdf

compiles the five papers from `paper-X/tex/` (three passes of `pdflatex`
each). A compiled `main.pdf` is committed in each `tex/` folder.

## Two tiers of computation

- **Quick checks** (`make check`): about 11 minutes on one core for the five
  papers together, about 1.1 GB of memory at the peak. They rerun the programs
  on the cases that finish in seconds or a few minutes and compare the output
  with the committed logs or with the expected lines stated in the READMEs.
- **Full runs** (`make full` prints them): the long computations behind the
  statements of the papers, from minutes to days and up to tens of gigabytes
  of memory. They are not rerun by `make check`; their commands, running
  times, memory use and committed logs are listed in each paper's `CLAIMS.md`
  and in the READMEs under `code/`.

## From results to claims

Each paper has a file `CLAIMS.md` with one row per computer-assisted
statement or confirmation mentioned in the paper's appendix: the statement
(numbered as in the compiled PDF), its status (proved by hand,
computer-assisted, or confirmation only), the program, the exact command, the
expected output or the committed log, the running time and memory, and the
tier.

- `paper-A/CLAIMS.md`
- `paper-B/CLAIMS.md`
- `paper-C/CLAIMS.md`
- `paper-D/CLAIMS.md`
- `paper-E/CLAIMS.md`

## Layout

    README.md, Makefile, CITATION.cff, LICENSE, LICENSE-data
    environment/        pinned versions, Dockerfile
    paper-A/ … paper-E/
        tex/            LaTeX sources, figures, compiled main.pdf
        code/           programs, READMEs and logs, as the paper's appendix describes them
        CLAIMS.md       statement → program → command → expected output
        Makefile        targets check, full, clean
    common/             a note: each paper's code folder is self-contained
    data/               a note on the large data files that are not in the repository
    .github/workflows/  runs `make check` on push and pull request

## How to cite

See `CITATION.cff`. A Zenodo DOI will be added at the first release. The
papers are in preparation; cite the paper that the code you use belongs to.

## Licences

- Code (all programs, Makefiles and scripts): MIT, see `LICENSE`.
- Logs and recorded outputs of the computations, data tables, and the paper
  sources and PDFs under `paper-X/tex/`: CC BY 4.0, see `LICENSE-data`.

## Use of AI tools

Large language models (Anthropic's Claude models) were used extensively in
this work, to develop proofs, to write and run the programs and their
independent re-implementations, to check proofs and computations, and to
draft the manuscripts; see the statement in each paper.
