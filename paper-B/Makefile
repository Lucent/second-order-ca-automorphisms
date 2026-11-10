# Paper B: Fronts at the speed of light.
#   make check   quick tier (about three minutes on one core); see CLAIMS.md
#   make full    prints the commands of the full tier and the logs they reproduce
#   make clean   removes the scratch directory .check/
#
# `make check` copies code/ to .check/code and runs everything there, so the
# committed outputs in code/ are never overwritten.  Outputs are compared with
# the committed files (diff) or with the expected lines stated in CLAIMS.md.

PYTHON ?= python3
CC ?= gcc
SHELL := /bin/sh
W := .check/code
REF := $(CURDIR)/code
comma := ,
FAIL = { echo "paper-B: MISMATCH in $(1)"; exit 1; }
HAVE_PYSAT := $(shell $(PYTHON) -c "import pysat" 2>/dev/null && echo yes)

.PHONY: check check-nosat check-sat full clean

check: check-nosat check-sat
ifeq ($(HAVE_PYSAT),yes)
	@echo "paper-B: all quick checks passed"
else
	@echo "paper-B: all quick checks that do not need PySAT passed (SAT checks skipped)"
endif

check-nosat:
	@rm -rf .check && mkdir .check && cp -r code .check/code
	@echo "paper-B [1/12] tp_direct.py: least P with (T_P) for all 65536 functions g at radius 2"
	@cd $(W) && rm -f out_leastP_r2.txt && $(PYTHON) tp_direct.py > new_tp_direct.txt \
	  && diff -q new_tp_direct.txt $(REF)/out_tp_direct.txt \
	  && diff -q out_leastP_r2.txt $(REF)/out_leastP_r2.txt || $(call FAIL,tp_direct.py)
	@echo "paper-B [2/12] scan_r2.py, scan_r2b.py: the second program for (T_P) at radius 2; Lambda_2, Lambda_3"
	@cd $(W) && $(PYTHON) scan_r2.py > new_scan_r2.txt \
	  && diff -q new_scan_r2.txt $(REF)/out_scan_r2.txt && cp new_scan_r2.txt out_scan_r2.txt \
	  && $(PYTHON) scan_r2b.py > new_scan_r2b.txt \
	  && diff -q new_scan_r2b.txt $(REF)/out_scan_r2b.txt || $(call FAIL,scan_r2.py / scan_r2b.py)
	@echo "paper-B [3/12] compare.py: the tables of least P of the two programs, for every g"
	@cd $(W) && rm -f out_leastP_r2_scan.txt && $(PYTHON) compare.py > new_compare.txt \
	  && diff -q new_compare.txt $(REF)/out_compare.txt \
	  && diff -q out_leastP_r2_scan.txt $(REF)/out_leastP_r2_scan.txt \
	  && grep -qx "lines: 65537 and 65537; differing lines: 0" new_compare.txt || $(call FAIL,compare.py)
	@echo "paper-B [4/12] xcheck.py: the two (T_P) routines on 153 functions; counts of (N), (D)"
	@cd $(W) && $(PYTHON) xcheck.py > new_xcheck.txt \
	  && diff -q new_xcheck.txt $(REF)/out_xcheck.txt || $(call FAIL,xcheck.py)
	@echo "paper-B [5/12] examples.py, checks3.py: the tables of examples at radius 1, 2, 3; Lambda_M"
	@cd $(W) && $(PYTHON) examples.py > new_examples.txt \
	  && diff -q new_examples.txt $(REF)/out_examples.txt || $(call FAIL,examples.py)
	@cd $(W) && $(PYTHON) checks3.py > new_checks3.txt \
	  && test "`sed -n '/== r=2 ==/,$$p' new_checks3.txt | grep -c 'Lambda_M SC M=0..3 \[True, True, True, True\]'`" = 8 \
	  || $(call FAIL,checks3.py)
	@cd $(W) && $(PYTHON) -c "import examples as e; print([all(e.lambda_strongly_connected(t,2,M) for M in range(4)) for t in e.RANDOM_R2], [all(e.lambda_strongly_connected(t,3,M) for M in range(3)) for t in e.RANDOM_R3])" > new_lambda.txt \
	  && grep -qx "\[True, True, True, True, True, True, True, True, True, True\] \[True, True, True, True, True, True\]" new_lambda.txt \
	  || $(call FAIL,Lambda_M for the rules of Tables 2 and 3)
	@echo "paper-B [6/12] ca.py: Lambda_M for the sixteen left-permutive elementary rules, M <= 5"
	@cd $(W) && $(PYTHON) ca.py > new_ca.txt \
	  && for r in 30 120 45 180 75 210 135 225; do \
	       grep -qx "$$r \[True, True, True, True, True, True\]" new_ca.txt || exit 1; done \
	  && grep -qx "90 \[True, False, False, False, False, False\]" new_ca.txt || $(call FAIL,ca.py)
	@echo "paper-B [7/12] lift_cert.py: the lift certificate for rules 30 and 45"
	@cd $(W) && $(PYTHON) lift_cert.py 30 2 1 50 > new_lift30.txt \
	  && grep -q "test(R_k) holds for all 0 <= k <= 2: True" new_lift30.txt \
	  && grep -q "repetition R_3 = R_2" new_lift30.txt || $(call FAIL,lift_cert.py 30 2 1 50)
	@cd $(W) && $(PYTHON) lift_cert.py 45 2 4 50 > new_lift45.txt \
	  && grep -q "test(R_k) holds for all 0 <= k <= 22: True" new_lift45.txt \
	  && grep -q "repetition R_23 = R_19" new_lift45.txt || $(call FAIL,lift_cert.py 45 2 4 50)
	@echo "paper-B [8/12] ref_check.py: the two-copy window graph (30 with w + M <= 6, except (5,1) and (6,0); 45, 75, 135 in eight cases)"
	@cd $(W) && ( for wm in "1 0" "1 1" "1 2" "1 3" "1 4" "1 5" "2 0" "2 1" "2 2" "2 3" "2 4" "3 0" "3 1" "3 2" "3 3" "4 0" "4 1" "4 2" "5 0"; do \
	    $(PYTHON) ref_check.py 30 $$wm > new_ref.txt || exit 1; \
	    grep -q "ok=True; .*strongly_conn=True; .*Hol(l_a,l_b) values=\[0, 1, 2, 3\]" new_ref.txt || exit 1; \
	  done ) || $(call FAIL,ref_check.py 30 W M)
	@cd $(W) && ( for r in 45 75 135; do for wm in "1 0" "2 0" "3 0" "4 0" "1 1" "2 1" "3 1" "1 2"; do \
	    $(PYTHON) ref_check.py $$r $$wm > new_ref.txt || exit 1; \
	    grep -q "ok=True; .*strongly_conn=True; .*Hol(l_a,l_b) values=\[0, 1, 2, 3\]" new_ref.txt || exit 1; \
	  done; done ) || $(call FAIL,ref_check.py RULE W M (two-copy window graph))
	@echo "paper-B [9/12] ref_check.py, ref_fast.py, ref_affine.py, ref_base.py: the holonomy behind the front"
	@cd $(W) && ( for r in 30 45 75 135 120 180 210 225; do for w in 2 3 4; do \
	    $(PYTHON) ref_check.py $$r $$w 0 > new_ref.txt || exit 1; \
	    grep -q "ok=True; .*strongly_conn=True; #Psi values of closed walks at a vertex = 32;" new_ref.txt || exit 1; \
	  done; done ) || $(call FAIL,ref_check.py RULE W 0)
	@cd $(W) && ( for wm in "5 0" "2 1" "3 1" "4 1" "2 2" "3 2" "2 3"; do \
	    $(PYTHON) ref_check.py 135 $$wm > new_ref.txt || exit 1; \
	    grep -q "ok=True; .*strongly_conn=True; #Psi values of closed walks at a vertex = 32;" new_ref.txt || exit 1; \
	  done ) || $(call FAIL,ref_check.py 135 W M)
	@cd $(W) && ( for wm in "6 0" "5 1" "4 2"; do \
	    $(PYTHON) ref_fast.py 135 $$wm > new_ref.txt || exit 1; \
	    grep -q "ok=True, 4|G|==|T|: True); .*SCCs=1; group generated by Psi on closed walks has order 32;" new_ref.txt || exit 1; \
	  done ) || $(call FAIL,ref_fast.py 135 W M)
	@cd $(W) && ( for r in 30 45 75 135 120 180 210 225; do \
	    $(PYTHON) ref_check.py $$r 3 0 flip > new_flip.txt || exit 1; \
	    test "`grep -c 'True$$' new_flip.txt`" = 3 && ! grep -q "False\|NOT" new_flip.txt || exit 1; \
	  done ) || $(call FAIL,ref_check.py RULE 3 0 flip)
	@cd $(W) && chk() { $(PYTHON) ref_affine.py $$1 $$2 0 > new_affine.txt || return 1; \
	      grep -q "SCCs=$$3, group generated by Psi over all closed walks: order $$4," new_affine.txt; }; \
	  chk 90 2 4 8 && chk 90 3 16 8 && chk 90 4 64 8 && chk 105 2 1 4 && chk 105 3 1 2 && chk 105 4 2 4 \
	  && chk 150 2 1 4 && chk 150 3 1 2 && chk 150 4 2 4 || $(call FAIL,ref_affine.py)
	@cd $(W) && $(PYTHON) ref_base.py base > new_base.txt \
	  && test "`grep -c 'closed 3-walk in G_2: True; Psi_2 = (T,l_a,l_b,N_a,N_b) = \[1, 0, 1, 1, 1\]' new_base.txt`" = 4 \
	  || $(call FAIL,ref_base.py base)
	@cd $(W) && $(PYTHON) ref_base.py lift 3 0 > new_lift.txt \
	  && grep -q "sampled: 61; lifted closed walk of the same length through zeta: 61" new_lift.txt || $(call FAIL,ref_base.py lift 3 0)
	@echo "paper-B [10/12] halfradius/hr_checks.py: the finite checks for the classes of 63 and 187"
	@cd $(W)/halfradius && $(PYTHON) hr_checks.py > new_hr_checks.txt \
	  && diff -q new_hr_checks.txt $(REF)/halfradius/out_hr_checks.txt || $(call FAIL,halfradius/hr_checks.py)

check-sat: check-nosat
ifeq ($(HAVE_PYSAT),yes)
	@echo "paper-B [11/12] sat/: commuting block maps (PySAT): rules 30, 45, 75, 135 at R <= 2; Table 1 at R = 1; strip"
	@cd $(W)/sat && $(PYTHON) centralizer.py --rules 30,45 --R 1,2 > new_cent.txt \
	  && test "`grep -c \"'other': 0\" new_cent.txt`" = 4 && test "`grep -c 'all named found: True' new_cent.txt`" = 4 \
	  && grep -q "rule 30R, radius R=1: .* 7 solutions" new_cent.txt && grep -q "rule 30R, radius R=2: .* 15 solutions" new_cent.txt \
	  && grep -q "rule 45R, radius R=1: .* 6 solutions" new_cent.txt && grep -q "rule 45R, radius R=2: .* 14 solutions" new_cent.txt \
	  || $(call FAIL,sat/centralizer.py --rules 30$(comma)45 --R 1$(comma)2)
	@cd $(W)/sat && chk() { $(PYTHON) cent_r.py $$1 $$2 $$3 > new_cent_r.txt || return 1; \
	      grep -q "\"nsol\": $$4, .*\"verified_fail\": 0, .*\"extra_automorphisms\": $$5, " new_cent_r.txt; }; \
	  chk ECA45 1 "" 6 0 && chk ECA135 1 "" 5 0 && chk ECA75 1 "" 6 0 && chk ECA45 2 "" 14 0 && chk ECA135 2 "" 13 0 \
	  && chk 0000000000110011 2 1 ".*" 0 && chk 0001000101110111 2 1 ".*" 3 \
	  && chk 0101010110101010 2 1 ".*" 0 && chk 0000111100001111 2 1 ".*" 5 || $(call FAIL,sat/cent_r.py)
	@cd $(W)/sat && $(PYTHON) invertible_sat.py --rules 45 --R 1 > new_inv.txt \
	  && grep -q "balanced commuting maps: 5;" new_inv.txt && grep -q "OTHER invertible: 0" new_inv.txt || $(call FAIL,sat/invertible_sat.py)
	@cd $(W)/sat && chk() { $(PYTHON) strip.py $$1 2 $$2 > new_strip.txt || return 1; grep -q "nonzero solutions: $$3 " new_strip.txt; }; \
	  chk 0000000000110011 1 0 && chk 0000000000110011 2 0 && chk 0001000101110111 1 1 && chk 0001000101110111 2 1 \
	  && chk 0101010110101010 1 15 && chk 0101010110101010 2 20 || $(call FAIL,sat/strip.py)
	@echo "paper-B [12/12] halfradius/kernel_sat.py: right-sided commuting maps (PySAT): 187, 119 (anticipation <= 6), 102"
	@cd $(W)/halfradius && ( for r in 187 119; do $(PYTHON) kernel_sat.py $$r 1,2,3,4,5,6 > new_kernel.txt || exit 1; \
	    test "`grep -c 'other solutions: 0 ' new_kernel.txt`" = 6 || exit 1; \
	    grep -q "anticipation <= 6: 13 solutions; powers F^b: b = \[-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6\]; other solutions: 0 " new_kernel.txt || exit 1; \
	  done ) || $(call FAIL,halfradius/kernel_sat.py)
	@cd $(W)/halfradius && $(PYTHON) kernel_sat.py 102 1,2,3 > new_kernel102.txt \
	  && grep -q "anticipation <= 1: 12 solutions; .* 6 not injective)" new_kernel102.txt \
	  && grep -q "anticipation <= 2: 48 solutions; .* 30 not injective)" new_kernel102.txt \
	  && grep -q "anticipation <= 3: 192 solutions; .* 114 not injective)" new_kernel102.txt || $(call FAIL,halfradius/kernel_sat.py 102)
else
	@echo "paper-B: PySAT (python-sat) not installed: skipping SAT checks [11/12] and [12/12] (sat/centralizer.py, sat/cent_r.py, sat/invertible_sat.py, sat/strip.py, halfradius/kernel_sat.py)"
endif

full:
	@echo "Full tier of paper B (not executed by make; one core; see CLAIMS.md):"
	@echo "  cd code"
	@echo "  python3 lift_cert.py 135 4 12 105        # 5 min, 1.2 GB; committed output: code/out_lift_cert_135.txt"
	@echo "  python3 ref_check.py 30 5 1              # 18 s, 2.5 GB"
	@echo "  python3 ref_check.py 30 6 0              # 81 s, 9.8 GB"
	@echo "  python3 ref_fast.py 135 6 1              # 18 s, 0.9 GB"
	@echo "  python3 ref_fast.py 135 7 0              # 35 s, 1.7 GB"
	@echo "  cd sat"
	@echo "  python3 centralizer.py --rules 30 --R 3  # 49 s, 1.2 GB: 27 solutions"
	@echo "  python3 cent_r.py 0000000000110011 2 2   # 22 s: nsol 11, extra_automorphisms 0"
	@echo "  python3 cent_r.py 0001000101110111 2 2   # 28 s: nsol 16, extra_automorphisms 7"
	@echo "  python3 cent_r.py 0101010110101010 2 2   # 87 s: nsol 64, extra_automorphisms 4"

clean:
	rm -rf .check
