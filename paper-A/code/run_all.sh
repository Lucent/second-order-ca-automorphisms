#!/bin/sh
# Runs the confirmations for the classes of 110, 14, 46, 24 and 189 and draws the figures; writes out_*.txt.
# The proofs for these classes use no computer; every run below is a confirmation.
set -e
cd "$(dirname "$0")"
gcc -O2 -Wall -o base_cases base_cases.c
t() { name=$1; shift; s=$(date +%s.%N); "$@" > out_$name.txt; e=$(date +%s.%N);
      printf '%s: %.1f s\n' "$name" "$(echo "$e - $s" | bc)"; tail -1 out_$name.txt; }
t base_cone  sh -c 'for r in 110 42 46; do for k in 1 2 3 4; do ./base_cases $r $k cone; done; done'
t base_short sh -c 'for k in 5 6 7 8 9 10; do ./base_cases 110 $k short; ./base_cases 42 $k short; done; for k in 4 5 6 7 8; do ./base_cases 46 $k short; done'
t checkA   python3 checkA_constructive.py 5 8 300
t checkB   python3 checkB_probes.py 5 1
t probes   python3 probes.py 1500
t checks14_46 python3 checks14_46.py
t checks24_rules python3 checks24_rules.py
t checks24_moves python3 checks24_moves.py
t checks24_ctrl  python3 checks24_ctrl.py
# figures (needs numpy and matplotlib): writes ../tex/figs/*.pdf (or ../figs/*.pdf) and out_figures.txt
if python3 -c 'import numpy, matplotlib' 2>/dev/null; then PY=python3
else PY="uv run --with numpy --with matplotlib python3"; fi
s=$(date +%s.%N); $PY figures.py > /dev/null; e=$(date +%s.%N)
printf 'figures: %.1f s\n' "$(echo "$e - $s" | bc)"; tail -1 out_figures.txt
