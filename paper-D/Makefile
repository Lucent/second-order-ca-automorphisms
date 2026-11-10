# Paper D: "The 88 second-order elementary cellular automata: a classification
# by automorphism group".
#
#   make check   quick tier: about one and a half minutes on one core.  Runs in
#                a scratch copy .check/code, so nothing in code/ is overwritten,
#                and compares every result with the committed logs.
#   make full    prints the commands of the long runs and the logs they reproduce
#   make clean   removes .check
#
# Requirements: Python 3 with NumPy and SciPy, a C compiler; the SAT steps need
# the package python-sat (PySAT) and are skipped, with a message, without it.

SHELL := /bin/bash
PYTHON ?= python3
CC ?= gcc
CFLAGS ?= -O2
export OMP_NUM_THREADS ?= 1
export OPENBLAS_NUM_THREADS ?= 1

# Absolute path of the interpreter when PYTHON is given as a relative path.
PY := $(if $(findstring /,$(PYTHON)),$(abspath $(PYTHON)),$(PYTHON))

HAVE_PYSAT := $(shell $(PY) -c "import pysat" 2>/dev/null && echo yes)

W := .check/code
# Replace running times (0.3s, [12s]) by a fixed token before comparing.
NT := sed -E 's/[0-9.]+s/Ts/g'
FAIL = { echo "paper-D: FAILED: $(1)"; exit 1; }

.PHONY: check full clean setup deps census census-lin census-sat searches windows reach

check: setup deps census census-lin census-sat searches windows reach
	@echo "paper-D: all quick checks passed"

setup:
	@rm -rf .check && mkdir .check && cp -r code .check/code

deps:
	@$(PY) -c "import numpy, scipy" 2>/dev/null \
	  || { echo "paper-D: NumPy and SciPy are required (pip install -r ../environment/requirements.txt)"; exit 1; }

# ---------------------------------------------------------------- census
census: setup
	@echo "paper-D census: classes88.py"
	@cd $(W)/census && $(PY) classes88.py > out_classes88.txt \
	  && diff out_classes88.txt logs/classes88.txt > /dev/null \
	  || $(call FAIL,classes88.py differs from logs/classes88.txt)

# The programs below import centralizer.py, which imports PySAT, although
# linear_90.py, rank90_packed.py, span90_R3.py and verify_90_structure.py do
# not call a solver.
ifeq ($(HAVE_PYSAT),yes)
census-lin: setup
	@echo "paper-D census: linear_90.py 1 (90R, radius 1: kernel dimension 5, 32 maps, 10 invertible)"
	@cd $(W)/census && $(PY) linear_90.py 1 > out_linear90.txt \
	  && $(NT) out_linear90.txt | diff - <(sed -n '1,5p' logs/linear_90_R12.txt | $(NT)) > /dev/null \
	  || $(call FAIL,linear_90.py 1 differs from the radius-1 block of logs/linear_90_R12.txt)
	@echo "paper-D census: rank90_packed.py 150 1|2|3 (kernel dimensions 4, 8, 12)"
	@cd $(W)/census && { $(PY) rank90_packed.py 150 1; $(PY) rank90_packed.py 150 2; $(PY) rank90_packed.py 150 3; } > out_rank150.txt \
	  && $(NT) out_rank150.txt | diff - <($(NT) logs/rank_150_R123.txt) > /dev/null \
	  || $(call FAIL,rank90_packed.py 150 differs from logs/rank_150_R123.txt)
	@echo "paper-D census: rank90_packed.py 90 3 (kernel dimension 65)"
	@cd $(W)/census && $(PY) rank90_packed.py 90 3 > out_rank90.txt \
	  && $(NT) out_rank90.txt | diff - <($(NT) logs/rank_90_R3.txt) > /dev/null \
	  || $(call FAIL,rank90_packed.py 90 3 differs from logs/rank_90_R3.txt)
	@echo "paper-D census: span90_R3.py, verify_90_structure.py"
	@cd $(W)/census && $(PY) span90_R3.py > out_span.txt \
	  && diff out_span.txt logs/span_90_R3.txt > /dev/null \
	  || $(call FAIL,span90_R3.py differs from logs/span_90_R3.txt)
	@cd $(W)/census && $(PY) verify_90_structure.py > out_verify90.txt \
	  && diff out_verify90.txt logs/verify_90_structure.txt > /dev/null \
	  || $(call FAIL,verify_90_structure.py differs from logs/verify_90_structure.txt)

census-sat: setup
	@set -e; cd $(W)/census; \
	  echo "paper-D census: centralizer.py --rules 110,30 --R 1,2 (7 and 15 solutions)"; \
	  $(PY) centralizer.py --rules 110,30 --R 1,2 > out_110_30.txt; \
	  $(NT) out_110_30.txt | grep -v '^$$' | sort > a.txt; \
	  awk '/radius R=3/{exit} {print}' logs/census_110_30_R123.txt | $(NT) | grep -v '^$$' | sort > b.txt; \
	  diff a.txt b.txt > /dev/null || $(call FAIL,centralizer.py 110/30 R=1/2 differs from logs/census_110_30_R123.txt); \
	  echo "paper-D census: centralizer.py --rules 54,90 --R 1 (13 and 32 solutions)"; \
	  $(PY) centralizer.py --rules 54,90 --R 1 > out_54_90.txt; \
	  $(NT) out_54_90.txt | sort > a.txt; $(NT) logs/census_54_90_R1.txt | sort > b.txt; \
	  diff a.txt b.txt > /dev/null || $(call FAIL,centralizer.py 54/90 R=1 differs from logs/census_54_90_R1.txt); \
	  echo "paper-D census: invertible_sat.py --rules 54 --R 2 --native (23 invertible maps)"; \
	  $(PY) invertible_sat.py --rules 54 --R 2 --native > out_inv54.txt; \
	  $(NT) out_inv54.txt | sort > a.txt; $(NT) logs/invertible_54_R2.txt | sort > b.txt; \
	  diff a.txt b.txt > /dev/null || $(call FAIL,invertible_sat.py 54 R=2 differs from logs/invertible_54_R2.txt); \
	  echo "paper-D census: crosscheck.py --R 1,2 --rules 110,30 --jobs 1 --rup (three solvers, proof check)"; \
	  $(PY) crosscheck.py --R 1,2 --rules 110,30 --jobs 1 --rup > out_crosscheck.txt; \
	  grep -q '^# 4 instances, 0 mismatches$$' out_crosscheck.txt || $(call FAIL,crosscheck.py reports mismatches); \
	  grep '^rule' out_crosscheck.txt | sed -E 's/ \[[0-9]+s\]$$//' | sort > a.txt; \
	  grep -E '^rule (110|30)R R=[12]:' logs/crosscheck_R123.txt | sed -E 's/ \[[0-9]+s\]$$//' | sort > b.txt; \
	  diff a.txt b.txt > /dev/null || $(call FAIL,crosscheck.py differs from logs/crosscheck_R123.txt)
else
census-lin census-sat: setup
	@echo "PySAT (python-sat) not installed: skipping SAT checks of paper-D ($@: every program of code/census except classes88.py)"
endif

# ---------------------------------------------------------------- searches
searches: setup
	@cd $(W)/searches && for p in w1 hb factor conserved; do $(CC) $(CFLAGS) -o $$p $$p.c || exit 1; done
	@echo "paper-D searches: w1 (rules 13, 110, 133)"
	@cd $(W)/searches && { ./w1 13 7 3 1; ./w1 110 9 2 1; ./w1 133 5 4 0; } | sed -E 's/ \(smallest of size [0-9]+\)//' > out_w1.txt \
	  && grep -E '^rule (13 k 7 p 3 s 1|110 k 9 p 2 s 1|133 k 5 p 4 s 0):' logs/w1_results.txt | diff - out_w1.txt > /dev/null \
	  || $(call FAIL,w1 differs from logs/w1_results.txt)
	@echo "paper-D searches: hb 110 10, hb 133 10"
	@cd $(W)/searches && { ./hb 110 10; ./hb 133 10; } > out_hb.txt \
	  && grep -E '^rule (110|133) k 10 ' logs/hb_results.txt | diff - out_hb.txt > /dev/null \
	  || $(call FAIL,hb differs from logs/hb_results.txt)
	@echo "paper-D searches: factor 46 3 3 k MODE, |k| <= 3"
	@cd $(W)/searches && for k in -3 -2 -1 0 1 2 3; do ./factor 46 3 3 $$k 1; ./factor 46 3 3 $$k 0; done > out_factor.txt \
	  && grep '^rule 46 r 3 p 3 ' logs/factor_results.txt | diff - out_factor.txt > /dev/null \
	  || $(call FAIL,factor differs from logs/factor_results.txt)
	@echo "paper-D searches: conserved RULE r, r = 1..4, rules 91 and 110; cq.py 91 2"
	@cd $(W)/searches && for f in 91 110; do for r in 1 2 3 4; do ./conserved $$f $$r; done; done > out_conserved.txt \
	  && { grep -E '^rule 91 r [1-4]:' logs/conserved_results.txt; grep -E '^rule 110 r [1-4]:' logs/conserved_results.txt; } | diff - out_conserved.txt > /dev/null \
	  || $(call FAIL,conserved differs from logs/conserved_results.txt)
	@cd $(W)/searches && $(PY) cq.py 91 2 | grep -q '^rule 91 r 2: solution dim 7, trivial rank 4, nontrivial 2$$' \
	  || $(call FAIL,cq.py 91 2)
	@echo "paper-D searches: nowalls110.py"
	@cd $(W)/searches && $(PY) nowalls110.py > out_nowalls110.txt \
	  && diff out_nowalls110.txt logs/nowalls110.txt > /dev/null \
	  || $(call FAIL,nowalls110.py differs from logs/nowalls110.txt)

# ---------------------------------------------------------------- windows
windows: setup
	@echo "paper-D windows: window_graphs.py 0,1,2,3,4 (all 88 classes, M <= 4)"
	@cd $(W)/windows && $(PY) window_graphs.py 0,1,2,3,4 > out_M0-4.txt \
	  && sed -E 's/([S.]{5})[S.]{2}/\1/g' window_graphs_M0-6.txt | diff - out_M0-4.txt > /dev/null \
	  || $(call FAIL,window_graphs.py differs from the columns M <= 4 of window_graphs_M0-6.txt)

# ---------------------------------------------------------------- reach
# moves_log.txt was written without the running time at the end of the lines
# "k ..." and with plain floats; both are removed before the comparison.
NM := sed -E 's/np\.float64\(([0-9.]+)\)/\1/g; s/ [0-9.]+s$$//'
reach: setup
	@echo "paper-D reach: reach1.py classes, reach1.py controls"
	@cd $(W)/reach && $(PY) reach1.py classes > out_classes.txt \
	  && diff out_classes.txt reach1_classes.txt > /dev/null \
	  || $(call FAIL,reach1.py classes differs from reach1_classes.txt)
	@cd $(W)/reach && $(PY) reach1.py controls > out_controls.txt \
	  && diff out_controls.txt reach1_controls.txt > /dev/null \
	  || $(call FAIL,reach1.py controls differs from reach1_controls.txt)
	@echo "paper-D reach: reach1_unique.py classes, reach1_unique.py controls"
	@cd $(W)/reach && $(PY) reach1_unique.py classes > out_uclasses.txt \
	  && diff out_uclasses.txt reach1_unique_classes.txt > /dev/null \
	  || $(call FAIL,reach1_unique.py classes differs from reach1_unique_classes.txt)
	@cd $(W)/reach && $(PY) reach1_unique.py controls > out_ucontrols.txt \
	  && diff out_ucontrols.txt reach1_unique_controls.txt > /dev/null \
	  || $(call FAIL,reach1_unique.py controls differs from reach1_unique_controls.txt)
	@echo "paper-D reach: moves.py for 114 (k <= 4), 159 (k <= 4), 110 (k <= 4)"
	@cd $(W)/reach && $(CC) $(CFLAGS) -shared -fPIC -o uf.so uf.c
	@cd $(W)/reach && $(PY) moves.py 114 --bx 1 --by 1 --speeds 1/1 --kmax 4 | $(NM) > out_114.txt \
	  && sed -n '236,257p' moves_log.txt | diff - out_114.txt > /dev/null \
	  || $(call FAIL,moves.py 114 differs from moves_log.txt)
	@cd $(W)/reach && $(PY) moves.py 159 --speeds 2/4 --kmax 4 | $(NM) > out_159.txt \
	  && sed -n '285,$$p' moves_log.txt | diff - out_159.txt > /dev/null \
	  || $(call FAIL,moves.py 159 differs from moves_log.txt)
	@cd $(W)/reach && $(PY) moves.py 110 --speeds 1/1 --kmax 4 | $(NM) > out_110.txt \
	  && sed -n '2,23p' moves_log.txt | diff - out_110.txt > /dev/null \
	  || $(call FAIL,moves.py 110 differs from the lines k <= 4 of moves_log.txt)

# ---------------------------------------------------------------- full
full:
	@echo "Long runs of paper D.  They are not executed here; each line gives the"
	@echo "directory, the command, the committed log it reproduces, and its cost."
	@echo
	@echo "code/census (needs python-sat):"
	@echo "  python3 centralizer.py --rules 110,30 --R 1,2,3            -> logs/census_110_30_R123.txt   (about 2 min)"
	@echo "  python3 centralizer.py --rules RULE --R 1,2,3 --max 5000   -> logs/census_R123_summary.txt  (72 rules; about 1 min per rule with few solutions, longer with many)"
	@echo "  python3 centralizer.py --rules 54 --R 2 --max 20000        -> logs/census_54_R2_full.txt    (cap reached after 127 s)"
	@echo "  python3 enum_r4.py RULE                                    -> logs/census_R4_summary.txt, logs/census_110_R4.txt, logs/census_30_R4.txt  (41 rules; 5-35 min and 5-7 GB per rule; writes raw_solutions_RULE_R4.txt)"
	@echo "  python3 identify_r4.py                                     -> logs/identify_R4.txt          (needs raw_solutions_{24,189,63,187}_R4.txt)"
	@echo "  python3 crosscheck.py --R 1,2,3 --jobs 12 --rup            -> logs/crosscheck_R123.txt      (200 instances; 4-8 min and 3 GB per instance at R = 3)"
	@echo "  python3 crosscheck.py --R 4 --rules 110,30 --solvers glucose4,cadical153 --jobs 2 -> logs/crosscheck_R4_110_30.txt (30 min, 6-15 GB per instance)"
	@echo "  python3 crosscheck.py --R 4 --solvers glucose4 --jobs 3    -> logs/crosscheck_R4.txt        (three invocations, recorded in the log; 10-30 min, 6-15 GB per instance)"
	@echo "  python3 linear_90.py 1,2                                   -> logs/linear_90_R12.txt"
	@echo "  python3 analyze_90.py 2 all                                -> logs/analyze_90_R2.txt        (5851 s)"
	@echo
	@echo "code/searches (38 rules each; see code/searches/README.md):"
	@echo "  ./w1 RULE K P S, (K,P) = (11,1),(9,2),(7,3),(5,4),(3,5), |S| <= P -> logs/w1_results.txt, logs/w1_smallest_sink.txt  (about a minute per line for K = 11)"
	@echo "  ./hb RULE 10                                               -> logs/hb_results.txt"
	@echo "  ./factor RULE r p k MODE, (r,p) = (5,1),(4,2),(3,3),(2,4),(1,5) -> logs/factor_results.txt"
	@echo "  ./conserved RULE r, r = 1..5                               -> logs/conserved_results.txt    (about a minute for r = 5)"
	@echo
	@echo "code/windows:"
	@echo "  python3 window_graphs.py                                   -> window_graphs_M0-6.txt        (20 min)"
	@echo
	@echo "code/reach:"
	@echo "  sh moves_all.sh                                            -> moves_log.txt                 (a few minutes for each of the four large runs)"

clean:
	rm -rf .check
