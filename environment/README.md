# Environment

## Versions

| component | version | used for |
|---|---|---|
| Python | 3.12 (CPython 3.12.14) | all Python programs |
| NumPy | 2.5.4 | most programs |
| python-sat (PySAT) | 1.9.dev21, with the bundled solvers CaDiCaL 1.5.3 (`cadical153`), MiniSat 2.2 (`minisat22`), Glucose 4 (`glucose4`) and MiniCard (`minicard`) | the SAT enumerations |
| SciPy | 1.18.1 | `paper-D/code/reach/reach1.py` only |
| Matplotlib | 3.11.2 | `paper-A/code/figures.py` only |
| gcc | 15.2.0 | the C programs (plain C99; `-O2`, and `-fopenmp` for the two OpenMP programs of the long runs) |

`requirements.txt` pins the four Python packages. The quick checks
(`make check` at the top level) were run with exactly these versions under
CPython 3.12.14 and gcc 15.2.0 on Linux (aarch64). The long runs whose logs
are committed were made with NumPy 2.5.4, python-sat 1.9.dev21 and gcc 15.2.0
on Linux (x86-64, 32 GB of memory); the statements of running time and memory
in each paper's `CLAIMS.md` and READMEs refer to that machine, one core unless
stated otherwise.

Programs that use neither a SAT solver nor NumPy need only the Python
standard library; each paper's `code/README.md` says what each program needs.

## Installing

    python3.12 -m venv .venv
    . .venv/bin/activate
    pip install -r environment/requirements.txt
    make check

## Docker

`Dockerfile` builds an image from `python:3.12-slim` with the pinned Python
packages and the C compiler of the base distribution (its gcc may be older
than 15.2; the C programs do not depend on the compiler version), and runs
`make check`:

    docker build -f environment/Dockerfile -t second-order-ca .
    docker run --rm second-order-ca
