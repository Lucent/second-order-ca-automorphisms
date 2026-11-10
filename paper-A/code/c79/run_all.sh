#!/bin/sh
# Runs the confirmations for the class of 79 (representative 93) and stores their output.
# The proof for this class uses no computer; every run below is a confirmation.
set -e
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
t() { name=$1; shift; s=$(date +%s); "$@" > out_$name.txt; e=$(date +%s);
      printf '%s (%d s): ' "$name" $((e-s)); tail -1 out_$name.txt; }
t check_F3_hand      python3 check_F3_hand.py
t check_F3_hand_core python3 check_F3_hand.py core
t c_F3               python3 c_F3.py
t check_F3_second    python3 check_F3_second.py
t check_H2H3         python3 check_H2H3.py 11 1500
t check_F4steps      python3 check_F4steps.py
t check_bisim_L3     python3 check_bisim.py 3
t check_bisim_L2     python3 check_bisim.py 2
t check_F4_graph     python3 check_F4_graph.py 3 4
t check_kinks        python3 check_kinks.py
t check_gluing       python3 check_gluing.py
t check_farglue      python3 check_farglue.py 1 40
t check_engine       python3 check_engine.py
t check_graphs_exact python3 check_graphs_exact.py
t check_closure64    python3 check_closure64.py
