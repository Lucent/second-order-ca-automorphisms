# Paper A: Reversible rule 110 has no hidden symmetries.
#   make check   quick tier (about a minute and a half on one core); see CLAIMS.md
#   make full    prints the commands of the full tier and the logs they reproduce
#   make clean   removes the scratch directory .check/
#
# `make check` copies code/ to .check/code and runs everything there, so the
# committed outputs in code/ are never overwritten.  Each output is compared
# with the committed file (lines "time ... s" excluded).

PYTHON ?= python3
CC ?= gcc
SHELL := /bin/sh
W := .check/code
REF := $(CURDIR)/code
FAIL = { echo "paper-A: MISMATCH in $(1)"; exit 1; }
# CMP new committed: equal up to the lines that report the running time
CMP = sh -c 'grep -v "^time " "$$0" > "$$0.a"; grep -v "^time " "$$1" > "$$0.b"; diff -q "$$0.a" "$$0.b" > /dev/null'
HAVE_NUMPY := $(shell $(PYTHON) -c "import numpy" 2>/dev/null && echo yes)
HAVE_MPL := $(shell $(PYTHON) -c "import numpy, matplotlib" 2>/dev/null && echo yes)

C79_QUICK := c_F3 check_F3_second check_kinks check_gluing check_engine check_graphs_exact check_closure64

.PHONY: check check-std check-c79 check-fig full clean

check: check-std check-c79 check-fig
ifeq ($(HAVE_MPL),yes)
	@echo "paper-A: all quick checks passed"
else
	@echo "paper-A: all quick checks that were not skipped passed (see the messages above)"
endif

check-std:
	@rm -rf .check && mkdir -p .check/tex/figs && cp -r code .check/code
	@echo "paper-A [1/8] base_cases.c, mode cone: orbits of the moves on the light-cone window, k <= 4 (110, 42, 46)"
	@$(CC) -O2 -Wall -o .check/base_cases code/base_cases.c
	@cd $(W) && ( for r in 110 42 46; do for k in 1 2 3 4; do ../base_cases $$r $$k cone || exit 1; done; done ) > new_base_cone.txt \
	  && $(CMP) new_base_cone.txt $(REF)/out_base_cone.txt || $(call FAIL,base_cases cone)
	@echo "paper-A [2/8] checkB_probes.py 5 1: every step of the proofs for 110 and 42, k <= 5"
	@cd $(W) && $(PYTHON) checkB_probes.py 5 1 > new_checkB.txt \
	  && $(CMP) new_checkB.txt $(REF)/out_checkB.txt || $(call FAIL,checkB_probes.py)
	@echo "paper-A [3/8] probes.py 1500: moves from probes and debris, bounce formula"
	@cd $(W) && $(PYTHON) probes.py 1500 > new_probes.txt \
	  && $(CMP) new_probes.txt $(REF)/out_probes.txt || $(call FAIL,probes.py)
	@echo "paper-A [4/8] checks14_46.py: conjugacy of 14 and 42, collision formula, bond bits of 46"
	@cd $(W) && $(PYTHON) checks14_46.py > new_checks14_46.txt \
	  && $(CMP) new_checks14_46.txt $(REF)/out_checks14_46.txt || $(call FAIL,checks14_46.py)
	@echo "paper-A [5/8] checks24_rules.py, checks24_moves.py, checks24_ctrl.py: classes of 24 and 189"
	@cd $(W) && ( for n in rules moves ctrl; do $(PYTHON) checks24_$$n.py > new_checks24_$$n.txt \
	  && $(CMP) new_checks24_$$n.txt $(REF)/out_checks24_$$n.txt || exit 1; done ) || $(call FAIL,checks24_*.py)

check-c79: check-std
ifeq ($(HAVE_NUMPY),yes)
	@echo "paper-A [6/8] c79/: check_F3_hand.py (twice), check_bisim.py 3 and 2, check_F4_graph.py 3 4"
	@cd $(W)/c79 && $(PYTHON) check_F3_hand.py > new_check_F3_hand.txt \
	  && $(CMP) new_check_F3_hand.txt $(REF)/c79/out_check_F3_hand.txt \
	  && $(PYTHON) check_F3_hand.py core > new_check_F3_hand_core.txt \
	  && $(CMP) new_check_F3_hand_core.txt $(REF)/c79/out_check_F3_hand_core.txt \
	  && $(PYTHON) check_bisim.py 3 > new_check_bisim_L3.txt \
	  && $(CMP) new_check_bisim_L3.txt $(REF)/c79/out_check_bisim_L3.txt \
	  && $(PYTHON) check_bisim.py 2 > new_check_bisim_L2.txt \
	  && $(CMP) new_check_bisim_L2.txt $(REF)/c79/out_check_bisim_L2.txt \
	  && $(PYTHON) check_F4_graph.py 3 4 > new_check_F4_graph.txt \
	  && $(CMP) new_check_F4_graph.txt $(REF)/c79/out_check_F4_graph.txt || $(call FAIL,c79 (F3 hand, bisim, F4 graph))
	@echo "paper-A [7/8] c79/: $(C79_QUICK)"
	@cd $(W)/c79 && ( for n in $(C79_QUICK); do $(PYTHON) $$n.py > new_$$n.txt \
	  && $(CMP) new_$$n.txt $(REF)/c79/out_$$n.txt || { echo "paper-A: MISMATCH in c79/$$n.py"; exit 1; }; done )
else
	@echo "paper-A: NumPy not installed: skipping the checks [6/8] and [7/8] for the class of 79 (code/c79)"
endif

check-fig: check-std
ifeq ($(HAVE_MPL),yes)
	@echo "paper-A [8/8] figures.py: the figures redrawn from exact simulations, with their checks"
	@cd $(W) && rm -f out_figures.txt && $(PYTHON) figures.py > /dev/null \
	  && $(CMP) out_figures.txt $(REF)/out_figures.txt \
	  && test `ls ../tex/figs/*.pdf | wc -l` = 6 || $(call FAIL,figures.py)
else
	@echo "paper-A: NumPy or matplotlib not installed: skipping [8/8] (figures.py)"
endif

full:
	@echo "Full tier of paper A (not executed by make; one core). The two scripts run every"
	@echo "confirmation, including the quick tier, and rewrite the out_*.txt files:"
	@echo "  cd code && ./run_all.sh        # about 5 min; the long steps are"
	@echo "      base_cases, mode short (70 s)        -> code/out_base_short.txt"
	@echo "      python3 checkA_constructive.py 5 8 300 (3 min) -> code/out_checkA.txt"
	@echo "  cd code/c79 && ./run_all.sh    # about 9 min; the long steps are"
	@echo "      python3 check_H2H3.py 11 1500 (3.7 min) -> code/c79/out_check_H2H3.txt"
	@echo "      python3 check_F4steps.py (3.5 min)      -> code/c79/out_check_F4steps.txt"
	@echo "      python3 check_farglue.py 1 40 (1.5 min) -> code/c79/out_check_farglue.txt"
	@echo "See CLAIMS.md."

clean:
	rm -rf .check
