#!/bin/sh
# Runs behind parts (2), (3) of the proposition on the reach of the edge-move method and the remark on
# the class of 31. Build the union-find helper first:  gcc -O2 -shared -fPIC -o uf.so uf.c
# Arguments of moves.py: RULE, the background (--bx, --by: its x-row and y-row over one spatial period;
# default the constant (0,0)), the speeds a/b of the kinks of Fix(sigma^a F^b) used as probes, and the
# range of k (number of cells next to the edge on which the moves are tabulated).
set -x
python moves.py 110 --speeds 1/1 --kmax 7
python moves.py 46  --speeds 1/1 --kmax 8
python moves.py 42  --speeds 1/1 --kmax 8
python moves.py 87  --speeds 1/3 --kmax 6
python moves.py 162 --bx 00 --by 01 --speeds 1/3 --kmax 4
python moves.py 25  --bx 01 --by 10 --speeds 1/1 --kmax 4
python moves.py 74  --speeds 1/1 --kmax 4
python moves.py 114 --bx 1 --by 1 --speeds 1/1 --kmax 4
python moves.py 191 --speeds 2/4 --kmax 4
python moves.py 159 --speeds 2/4 --kmax 4
