# common

There is no shared code in this directory, by design.

Each paper's appendix names the programs of its own `code/` folder, and each
`paper-X/code/` is self-contained: it can be copied out of the repository and
run as its README describes. Some files therefore occur in more than one
folder (for example `centralizer.py`, the SAT enumeration of commuting block
maps, in `paper-B/code/sat/`, `paper-C/code/`, `paper-D/code/census/` and
`paper-E/code/census/`). The copies are kept as each paper describes them and
are not required to be identical.
