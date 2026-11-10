# Paper E: "The cellular automata that commute with Rule 54, and rigidity of
# the second-order automaton 94R".
#
#   make check   quick tier: a few minutes on one core.  Runs in a scratch copy
#                .check/code, so nothing in code/ is overwritten, and compares
#                every result with the expected output stated in the READMEs of
#                code/ or with the committed logs.
#   make full    prints the commands of the long runs and the logs they reproduce
#   make clean   removes .check
#
# Requirements: Python 3 with NumPy, a C compiler with OpenMP (gcc); the SAT
# steps need the package python-sat (PySAT) and are skipped, with a message,
# without it.

SHELL := /bin/bash
PYTHON ?= python3
CC ?= gcc
export OMP_NUM_THREADS ?= 1
export OPENBLAS_NUM_THREADS ?= 1

# Absolute path of the interpreter when PYTHON is given as a relative path.
PY := $(if $(findstring /,$(PYTHON)),$(abspath $(PYTHON)),$(PYTHON))
HAVE_PYSAT := $(shell $(PY) -c "import pysat" 2>/dev/null && echo yes)

W := .check/code
# Replace running times (0.3s, [12s]) by a fixed token before comparing.
NT := sed -E 's/[0-9.]+s/Ts/g'
FAIL = { echo "paper-E: FAILED: $(1)"; exit 1; }
# Shell functions used in the recipes: "has FILE LINE" (the file contains
# exactly this line) and "sub FILE TEXT" (the file contains this text).
H := has(){ grep -qxF -- "$$2" "$$1" || { echo "paper-E: FAILED: expected line not found in $$1: $$2"; exit 1; }; }; \
     sub(){ grep -qF -- "$$2" "$$1" || { echo "paper-E: FAILED: expected text not found in $$1: $$2"; exit 1; }; };

.PHONY: check full clean setup deps closure kink lemmaL rule155 rule91 rule94 rule95 rule95-sat census sh translation95

check: setup deps closure kink lemmaL rule155 rule91 rule94 rule95 rule95-sat census sh translation95
	@echo "paper-E: all quick checks passed"

setup:
	@rm -rf .check && mkdir .check && cp -r code .check/code

deps:
	@$(PY) -c "import numpy" 2>/dev/null \
	  || { echo "paper-E: NumPy is required (pip install -r ../environment/requirements.txt)"; exit 1; }

# ---------------------------------------------------------------- closure
closure: setup
	@echo "paper-E closure: dc RULE R for 91, 94, 155 and R = 1, 2, 3"
	@$(H) cd $(W)/closure && $(CC) -O3 -o dc dc.c \
	  && for f in 91 94 155; do for R in 1 2 3; do ./dc $$f $$R | grep -E '^(rule|COMPLETE|STALLED)'; done; done > out_dc.txt \
	  && awk '/ R 4 /{skip=1; next} skip{skip=0; next} {print}' logs/dc_R1-4.log | diff - out_dc.txt > /dev/null \
	  || $(call FAIL,dc differs from the lines R <= 3 of logs/dc_R1-4.log)

# ---------------------------------------------------------------- kink
kink: setup
	@echo "paper-E kink: t_pair.py, t_R1.py 8, t_dal.py, t_train.py, t_step4.py, t95a.py, t95b.py 9"
	@$(H) cd $(W)/kink && $(PY) t_pair.py > out_t_pair.txt && test `wc -l < out_t_pair.txt` -eq 8 \
	  || $(call FAIL,t_pair.py)
	@$(H) cd $(W)/kink && $(PY) t_R1.py 8 > out_t_R1.txt \
	  && sed -E 's/cases (266000|466930)/cases N/' logs/t_R1_10.log | diff - <(sed -E 's/cases (28080|49192)/cases N/' out_t_R1.txt) > /dev/null \
	  && test `grep -c 'cases 28080 count mismatches 0 location mismatches 0' out_t_R1.txt` -eq 2 \
	  && test `grep -c 'cases 49192 count mismatches 0 location mismatches 0' out_t_R1.txt` -eq 2 \
	  || $(call FAIL,t_R1.py 8)
	@$(H) cd $(W)/kink && $(PY) t_dal.py > out_t_dal.txt && has out_t_dal.txt 'cases 800 bad 0'
	@$(H) cd $(W)/kink && $(PY) t_train.py > out_t_train.txt \
	  && grep -q '^beta 0 .*ok True' out_t_train.txt && ! grep '^beta 0 ' out_t_train.txt | grep -q 'ok False' \
	  || $(call FAIL,t_train.py)
	@$(H) cd $(W)/kink && $(PY) t_step4.py > out_t_step4.txt \
	  && sub out_t_step4.txt 'cores 351' && sub out_t_step4.txt 'D0 cores passing both shifts: [()] 1'
	@$(H) cd $(W)/kink && $(PY) t95a.py > out_t95a.txt && test `wc -l < out_t95a.txt` -eq 4 \
	  || $(call FAIL,t95a.py)
	@$(H) cd $(W)/kink && $(PY) t95b.py 9 > out_t95b.txt \
	  && sub out_t95b.txt 'cores 40' \
	  && sub out_t95b.txt 'first-row tests 6000 mismatch 0 second-row tests 6000 mismatch 0' \
	  && sub out_t95b.txt 'solutions of both relations: [((0,), (0,), 11, 10)]'

# ---------------------------------------------------------------- lemmaL
lemmaL: setup
	@echo "paper-E lemmaL: climb 22 91|94, climb_gains, dalembert.py, separ.py"
	@$(H) cd $(W)/lemmaL && $(CC) -O2 -o climb climb.c && $(CC) -O2 -o climb_gains climb_gains.c
	@$(H) cd $(W)/lemmaL && { ./climb 22 91; ./climb 22 94; } > out_climb.txt \
	  && has out_climb.txt 'rule 91 L=22 steps 8388608  viol(i)=0 viol(iii)=0 viol(ii)=0' \
	  && has out_climb.txt 'rule 94 L=22 steps 8388608  viol(i)=0 viol(iii)=0 viol(ii)=0'
	@$(H) cd $(W)/lemmaL && for r in 91 94 95; do ./climb_gains $$r 20 6; done > out_gains.txt \
	  && test "`grep -o 'min 2-step gain [0-9]*' out_gains.txt | tr '\n' ' '`" = "min 2-step gain 4 min 2-step gain 3 min 2-step gain 4 " \
	  && test `grep -c 'monotone violations 0' out_gains.txt` -eq 3 \
	  || $(call FAIL,climb_gains)
	@$(H) cd $(W)/lemmaL && $(PY) dalembert.py > out_dalembert.txt \
	  && for r in 91 94 95 155; do grep -q "^$$r formula violations: (0, 14400)" out_dalembert.txt || $(call FAIL,dalembert.py rule $$r); done
	@$(H) cd $(W)/lemmaL && { $(PY) separ.py 91 7 40 4; $(PY) separ.py 94 7 40 4; } > out_separ.txt \
	  && test `grep -c 'separated by T=40: 16384/16384' out_separ.txt` -eq 5 \
	  && test `grep -c 'left S_R after entering: 0' out_separ.txt` -eq 5 \
	  && test `grep -c 'firings at t>=T/2: 0' out_separ.txt` -eq 5 \
	  || $(call FAIL,separ.py)

# ---------------------------------------------------------------- rule155
rule155: setup
	@echo "paper-E rule155: edge155.py, thmW.py, checkW.py, stripsim.py, strip.py, s4.py"
	@$(H) cd $(W)/rule155 && $(PY) edge155.py > out_edge.txt \
	  && has out_edge.txt 'trials 3000 violations 0 11-front cases firing every step 799'
	@$(H) cd $(W)/rule155 && $(PY) thmW.py 1 3000 > out_thmW.txt \
	  && has out_thmW.txt "checks 99059 violations {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'phi': 0}"
	@$(H) cd $(W)/rule155 && $(PY) checkW.py 40 4 3 > out_checkW.txt \
	  && has out_checkW.txt 'W-check T 40 wmax 4 violations 0 checks 40920'
	@$(H) cd $(W)/rule155 && $(PY) stripsim.py > out_stripsim.txt \
	  && has out_stripsim.txt 'recursion violations 0 genuine fronts 150 late strip not on a tree cycle 0'
	@$(H) cd $(W)/rule155 && $(PY) strip.py 1 40 > out_strip.txt \
	  && has out_strip.txt "line 1 branches 1 {'double': 1} periods [2]" \
	  && has out_strip.txt "line 7 branches 2 {'free': 1} periods [4]" \
	  && has out_strip.txt "line 11 branches 2 {'slaved': 1, 'double': 1} periods [4, 8]" \
	  && has out_strip.txt "line 39 branches 34 {'slaved': 29, 'double': 1, 'free': 2} periods [8, 16]"
	@$(H) cd $(W)/rule155 && $(PY) s4.py 18 > out_s4.txt \
	  && sub out_s4.txt 'width 18 00-front cycles 220' \
	  && sub out_s4.txt 'Phi- words 86 shifted 11-fronts 106 shifted 11-fronts with the final 11 on the last two lines 28 other 0'
	@echo "paper-E rule155: mkbase.py and lc (window closure, R = 2, 3)"
	@$(H) cd $(W)/rule155 && $(CC) -O3 -fopenmp -o lc lc.c \
	  && $(PY) mkbase.py 155 2 "1,1;-1,1" base2.bin > out_lc.txt && ./lc dc 155 2 2 base2.bin >> out_lc.txt 2>&1 \
	  && $(PY) mkbase.py 155 3 "1,1;-1,1" base3.bin >> out_lc.txt && ./lc dc 155 3 2 base3.bin >> out_lc.txt 2>&1 \
	  && sub out_lc.txt 'unforced 0 iterations 10' && sub out_lc.txt 'unforced 0 iterations 18'
	@echo "paper-E rule155: sonly.py (R = 2, 3), mkseed.py and sc (R = 4)"
	@$(H) cd $(W)/rule155 && { $(PY) sonly.py 2 none; $(PY) sonly.py 3 none; $(PY) sonly.py 2 fronts 01 1 -2 2; \
	    $(PY) sonly.py 3 fronts 01 1 -3 3; $(PY) sonly.py 2 fronts 01 1 -2 1; } > out_sonly.txt \
	  && sub out_sonly.txt 'unforced (968, 7)' && sub out_sonly.txt 'unforced (16140, 10)' \
	  && sub out_sonly.txt 'unforced (0, 15)' && sub out_sonly.txt 'unforced (0, 31)' \
	  && sub out_sonly.txt 'unforced (875,'
	@$(H) cd $(W)/rule155 && $(CC) -O3 -o sc sc.c \
	  && { $(PY) mkseed.py 4 01 1 -4 5 0 s4.bin; ./sc 4 s4.bin; $(PY) mkseed.py 4 - 1 0 0 0 s4n.bin; ./sc 4 s4n.bin; } > out_sc.txt \
	  && sub out_sc.txt 'RESULT R 4 seed 276 unforced 0' && sub out_sc.txt 'RESULT R 4 seed 54 unforced 261168'

# ---------------------------------------------------------------- rule91
rule91: setup
	@echo "paper-E rule91: t91b.py 12, t91c.py 12, s91.py, cq_check.py, q1_identity.py, cq.py 91 2, q91.py"
	@$(H) cd $(W)/rule91 && $(PY) t91b.py 12 > out_t91b.txt \
	  && sub out_t91b.txt '(1) cores of Phi_91 of length <= 12 : 778  transparent to the bare defect: 41  disagreements with 3-sparse: 0 []' \
	  && sub out_t91b.txt '(2) first-ray lemma: cases 4000 violations 0 ; second-ray parity lemma: cases 1995 violations 0'
	@$(H) cd $(W)/rule91 && $(PY) t91c.py 12 > out_t91c.txt \
	  && sub out_t91c.txt 'cores of length <= 12 : 778 transparent: 41 disagreements with 3-sparse: 0'
	@$(H) cd $(W)/rule91 && $(PY) s91.py > out_s91.txt \
	  && sub out_s91.txt 'shape changed ((0,), (0, 0), 0, (4, (0, 0)), (4, (0, 1, 0)), 4, 4)'
	@$(H) cd $(W)/rule91 && for r in 91 94 95; do $(PY) cq_check.py $$r; done > out_cq_check.txt \
	  && sub out_cq_check.txt 'rule 91 Q1 conserved: True  Q2 conserved: True' \
	  && sub out_cq_check.txt 'rule 94 Q1 conserved: False  Q2 conserved: True' \
	  && sub out_cq_check.txt 'rule 95 Q1 conserved: False  Q2 conserved: True'
	@$(H) cd $(W)/rule91 && $(PY) q1_identity.py > out_q1.txt \
	  && sub out_q1.txt 'rule 91 regroup violations 0 description of s violations 0 termwise invariance violations 0 | periodic configs 87380 Q1 violations 0'
	@$(H) cd $(W)/rule91 && $(PY) cq.py 91 2 > out_cq.txt \
	  && sub out_cq.txt 'nullity 7' && sub out_cq.txt 'rank trivial 4 rank all h 6 nontrivial 2'
	@$(H) cd $(W)/rule91 && for r in 91 94 95; do $(PY) q91.py $$r; done > out_q91.txt \
	  && sub out_q91.txt 'rule 91 configs 24844 Q1 violations 0 Q2 violations 0' \
	  && sub out_q91.txt 'rule 94 configs 24844 Q1 violations 15580 Q2 violations 0' \
	  && sub out_q91.txt 'rule 95 configs 24844 Q1 violations 13296 Q2 violations 0'

# ---------------------------------------------------------------- rule94
rule94: setup
	@echo "paper-E rule94: verify1.py (lib 300, direct 200), verify2.py 6 8, verify3.py, closure_bf.py 4 9 9 2"
	@$(H) cd $(W)/rule94 && $(PY) verify1.py lib 300 > out_verify1.txt && $(PY) verify1.py direct 200 >> out_verify1.txt \
	  && test `grep -c ' bad 0$$' out_verify1.txt` -eq 14 && ! grep ' bad ' out_verify1.txt | grep -qv ' bad 0$$' \
	  && has out_verify1.txt 'seed 28 bad 0' \
	  || $(call FAIL,verify1.py)
	@$(H) cd $(W)/rule94 && $(PY) verify2.py 6 8 > out_verify2.txt \
	  && has out_verify2.txt "core length 8 objects certified so far 86 {'passes': 451, 'maxpass': 8}"
	@$(H) cd $(W)/rule94 && $(PY) verify3.py > out_verify3.txt \
	  && has out_verify3.txt 'points checked 51200 bad 0' && has out_verify3.txt 'part (ii): points checked 20800 bad 0'
	@$(H) cd $(W)/rule94 && $(PY) closure_bf.py 4 9 9 2 > out_closure_bf.txt \
	  && has out_closure_bf.txt 'round 5 known 153 of 153 ; core<= 9 : 153 / 153' && has out_closure_bf.txt 'missing []'
	@echo "paper-E rule94: ta.py 94 13, tc.py (10 12; 4 8), te.py 94 8 1,2,3,4,5,6, te.py 94 9 4,6,8 123"
	@$(H) cd $(W)/rule94 && $(PY) ta.py 94 13 > out_ta.txt \
	  && has out_ta.txt 'rule 94 objects 2513' && test `grep -c 'fail=0 *$$' out_ta.txt` -eq 39 && test `wc -l < out_ta.txt` -eq 40 \
	  && grep -qE '^lem B output +n=810 +fail=0 *$$' out_ta.txt \
	  && grep -qE '^prepend c=0: S\(.* n=10810 +fail=0 *$$' out_ta.txt && grep -qE '^prepend c=1: S\(.* n=14320 +fail=0 *$$' out_ta.txt \
	  && grep -qE '^prependZ u=d\+3: S=.* n=2864 +fail=0 *$$' out_ta.txt && grep -qE '^prependZ u>=d\+4: S=.* n=11456 +fail=0 *$$' out_ta.txt \
	  || $(call FAIL,ta.py 94 13)
	@$(H) cd $(W)/rule94 && { $(PY) tc.py 10 12; $(PY) tc.py 4 8; } > out_tc.txt \
	  && has out_tc.txt "rule 94 L 10 maxcore 12 objects 1434 {'bare': 4, 'prepend': 1079, 'prependZ': 350, 'zero': 1} scatterings 58164 largest core seen 45 FAILURES 0" \
	  && has out_tc.txt "rule 94 L 4 maxcore 8 objects 153 {'bare': 4, 'prepend': 112, 'prependZ': 36, 'zero': 1} scatterings 2474 largest core seen 22 FAILURES 0"
	@$(H) cd $(W)/rule94 && $(PY) te.py 94 8 1,2,3,4,5,6 > out_te.txt \
	  && sub out_te.txt 'scattering table: 23850 simulations' && test `grep -c 'closure has 153 of 153' out_te.txt` -eq 6 \
	  || $(call FAIL,te.py 94 8)
	@$(H) cd $(W)/rule94 && $(PY) te.py 94 9 4,6,8 123 > out_te9.txt \
	  && test `grep -c 'closure has 4 of 267' out_te9.txt` -eq 3 || $(call FAIL,te.py 94 9 with (C0) left out)

# ---------------------------------------------------------------- rule95
rule95: setup
	@echo "paper-E rule95: core.py, gen7.py, verify7 (reduced), props7.py, standalone7.py"
	@$(H) cd $(W)/rule95 && $(PY) core.py > out_core.txt && test `grep -c 'mismatches 0$$' out_core.txt` -eq 3 \
	  || $(call FAIL,core.py)
	@$(H) cd $(W)/rule95 && mkdir gen && $(PY) gen7.py gen > out_gen7.txt \
	  && has out_gen7.txt 'saved; ones: 16384 16384' \
	  && cmp -s gen/tx7.bin tx7.bin && cmp -s gen/ty7.bin ty7.bin && cmp -s gen/tx7.npy tx7.npy && cmp -s gen/ty7.npy ty7.npy \
	  || $(call FAIL,gen7.py does not reproduce the supplied tables)
	@$(H) cd $(W)/rule95 && echo "5a8b4dc6a76c188f2d5ccab5b224ddc3  tx7.bin" | md5sum -c --quiet - \
	  && echo "fda0396a76179ee07663e93e84d27166  ty7.bin" | md5sum -c --quiet - \
	  || $(call FAIL,md5 sums of tx7.bin and ty7.bin)
	@$(H) cd $(W)/rule95 && $(CC) -O2 -fopenmp -o verify7 verify7.c \
	  && ./verify7 tx7.bin ty7.bin -1 22 22 > out_verify7.txt \
	  && has out_verify7.txt '(A) class 0 (y-bit): 2^17 windows, mismatches 0' \
	  && has out_verify7.txt '(A) class 1 (x-bit): 2^17 windows, mismatches 0' \
	  && has out_verify7.txt '(B) class 0 (y-bit): 2^22 windows, violations of (sigma h)^2=id: 0' \
	  && has out_verify7.txt '(B) class 1 (x-bit): 2^22 windows, violations of (sigma h)^2=id: 0' \
	  && has out_verify7.txt '(C) 2^22 fully random windows: commutation mismatches 0, involution violations 0'
	@$(H) cd $(W)/rule95 && ./verify7 tx7.bin ty7.bin 777 22 22 > out_verify7neg.txt \
	  && has out_verify7neg.txt '(A) class 0 (y-bit): 2^17 windows, mismatches 8' \
	  && has out_verify7neg.txt '(A) class 1 (x-bit): 2^17 windows, mismatches 4' \
	  && has out_verify7neg.txt '(B) class 0 (y-bit): 2^22 windows, violations of (sigma h)^2=id: 256' \
	  && has out_verify7neg.txt '(B) class 1 (x-bit): 2^22 windows, violations of (sigma h)^2=id: 128'
	@$(H) cd $(W)/rule95 && $(PY) props7.py > out_props7.txt \
	  && sub out_props7.txt 'restriction to Phi^+ = sigma^-1 o (phi_* on each class): 3025 words, 0 failures; restriction to Phi^- = sigma^-1: 0 failures' \
	  && sub out_props7.txt 'period 5: 1024 configurations, failures of hF=Fh or (sigma h)^2=id: 0; sigma h moves 310; (n<=7) h(c) outside the <sigma,F>-orbit of c: 310' \
	  && sub out_props7.txt 'period 7: 16384 configurations, failures of hF=Fh or (sigma h)^2=id: 0; sigma h moves 7168; (n<=7) h(c) outside the <sigma,F>-orbit of c: 7168' \
	  && sub out_props7.txt 'period 9: 262144 configurations, failures of hF=Fh or (sigma h)^2=id: 0; sigma h moves 144396' \
	  && test `grep -c 'failures of hF=Fh or (sigma h)^2=id: 0;' out_props7.txt` -eq 9 \
	  || $(call FAIL,props7.py)
	@$(H) cd $(W)/rule95 && $(PY) standalone7.py tx7.npy ty7.npy > out_standalone7.txt \
	  && sub out_standalone7.txt 'F-period of c: 15' \
	  && sub out_standalone7.txt 'h(F(d)) == F(h(d))      on the F-orbit: True' \
	  && sub out_standalone7.txt '(sigma h)^2 (d) == d    on the F-orbit: True' \
	  && sub out_standalone7.txt 'sigma h (d) != d        on the F-orbit: True' \
	  && sub out_standalone7.txt 'h(c) = sigma^a F^b (c) for some a, b  : False'
	@echo "paper-E rule95: enum_golden.py 7"
	@$(H) cd $(W)/rule95 && $(PY) enum_golden.py 7 > out_enum_golden.txt \
	  && head -n 6 out_enum_golden.txt | diff - <(head -n 6 logs/enum_golden_8.log) > /dev/null \
	  && test `wc -l < out_enum_golden.txt` -eq 6 \
	  || $(call FAIL,enum_golden.py 7 differs from the lines w <= 7 of logs/enum_golden_8.log)
	@echo "paper-E rule95: census_cmp.py 1|2|3, block_test.py 4, window_graph.py"
	@$(H) cd $(W)/rule95 && { $(PY) census_cmp.py 1; $(PY) census_cmp.py 2; } > out_census_cmp12.txt \
	  && has out_census_cmp12.txt "R=1: {'total': 6, 'commute': 6, 'const': 1, 'sublattice': 5, 'form': 5, 'conserving': 5, 'both_injective': 5, 'noninjective': 0, 'table_equal': 5, 'P': 1, 'S': 4} distinct (phi,psi,n): 5" \
	  && has out_census_cmp12.txt "R=2: {'total': 14, 'commute': 14, 'const': 1, 'sublattice': 13, 'form': 13, 'conserving': 13, 'both_injective': 13, 'noninjective': 0, 'table_equal': 13, 'P': 9, 'S': 4} distinct (phi,psi,n): 13"
	@$(H) cd $(W)/rule95 && $(PY) census_cmp.py 3 > out_census_cmp3.txt \
	  && diff out_census_cmp3.txt logs/census_cmp_R3.log > /dev/null \
	  || $(call FAIL,census_cmp.py 3 differs from logs/census_cmp_R3.log)
	@$(H) cd $(W)/rule95 && $(PY) block_test.py 4 > out_block.txt \
	  && has out_block.txt 'length 4: 256 words x 4 positions; defect witness missing for 0 (word,position) pairs, largest first time of difference 35; random pair of contexts: missing 0, largest first time 15'
	@$(H) cd $(W)/rule95 && $(PY) window_graph.py > out_wg.txt \
	  && has out_wg.txt 'rule 95R, W^(+1,1)_M strongly connected for M=0..7: [True, True, True, True, True, True, True, True]' \
	  && has out_wg.txt 'rule 95R, W^(-1,1)_M strongly connected for M=0..7: [True, True, True, True, True, True, True, True]' \
	  && has out_wg.txt 'rule 95R, W^(+0,1)_M strongly connected for M=0..7: [True, True, True, True, True, True, True, True]' \
	  && has out_wg.txt 'rule 90R, W^(+1,1)_M strongly connected for M=0..4: [True, False, False, False, False]' \
	  && has out_wg.txt 'rule 90R, W^(-1,1)_M strongly connected for M=0..4: [True, False, False, False, False]' \
	  && has out_wg.txt 'rule 90R, W^(+0,1)_M strongly connected for M=0..4: [True, True, True, True, True]'
	@echo "paper-E rule95: verify36.py (the 36 automorphisms of radius 7)"
	@$(H) cd $(W)/rule95 && $(PY) verify36.py > out_verify36.txt \
	  && has out_verify36.txt 'tables verified: 36 all ok: True ; pairwise distinct: True' \
	  && has out_verify36.txt 'SAT search: 36 automorphisms of sublattice type S and radius <= 7 outside <sigma,F>; equal to the set of the 36 tables above: True' \
	  && test `grep -c 'commutation mismatches 0, restrictions as predicted: True, neighbourhood \[-7,7\]' out_verify36.txt` -eq 36 \
	  || $(call FAIL,verify36.py)

ifeq ($(HAVE_PYSAT),yes)
rule95-sat: setup rule95
	@echo "paper-E rule95 (SAT): sat_sublattice.py R all and R aut (R <= 5), with CaDiCaL"
	@$(H) cd $(W)/rule95 && for R in 1 2 3 4 5; do $(PY) sat_sublattice.py $$R all; done > out_sat_all.txt \
	  && diff <($(NT) out_sat_all.txt) <($(NT) logs/sat_sublattice_all_R1-5.log) > /dev/null \
	  || $(call FAIL,sat_sublattice.py R all differs from logs/sat_sublattice_all_R1-5.log)
	@$(H) cd $(W)/rule95 && for R in 1 2 3 4 5; do $(PY) sat_sublattice.py $$R aut; done > out_sat_aut.txt \
	  && diff <($(NT) out_sat_aut.txt) <(awk '/^R=6 /{exit} {print}' logs/sat_sublattice_aut_R1-6.log | $(NT)) > /dev/null \
	  || $(call FAIL,sat_sublattice.py R aut differs from the lines R <= 5 of logs/sat_sublattice_aut_R1-6.log)
	@echo "paper-E rule95 (SAT): sat_post.py 4 all, sat_post.py 5 aut, verify_sub.py 3"
	@$(H) cd $(W)/rule95 && $(PY) sat_post.py 4 all > out_sat_post.txt && $(PY) sat_post.py 5 aut >> out_sat_post.txt \
	  && sub out_sat_post.txt 'R=4 type P: 226 solutions; distinct class maps occurring in restrictions: 16, conserving 15, zero map 1' \
	  && sub out_sat_post.txt 'R=4 type S: 37 solutions; distinct class maps occurring in restrictions: 10, conserving 9, zero map 1' \
	  && sub out_sat_post.txt 'found 25, named 25, others 0' && sub out_sat_post.txt 'found 36, named 36, others 0'
	@$(H) cd $(W)/rule95 && $(PY) verify_sub.py 3 sat_R3_all.pkl > out_verify_sub.txt \
	  && has out_verify_sub.txt "R=3 type P: {'n': 9, 'const': 1, 'form': 9, 'conserving': 9, 'inj': 9, 'table_equal': 9}" \
	  && has out_verify_sub.txt "R=3 type S: {'n': 36, 'const': 1, 'form': 36, 'conserving': 36, 'inj': 16, 'table_equal': 36}"
	@echo "paper-E rule95 (SAT): rule 175, R <= 5, two encoders; rule 250, R = 3"
	@$(H) cd $(W)/rule95 && for R in 1 2 3 4 5; do $(PY) sat_sublattice.py $$R all cadical153 175; done > out_175.txt \
	  && diff <($(NT) out_175.txt) <(head -n `wc -l < out_175.txt` logs/sat_sublattice_rule175.log | $(NT)) > /dev/null \
	  && test `grep -c "both types together {'named': [0-9]*, 'const': 0, 'other': 0}" out_175.txt` -eq 5 \
	  || $(call FAIL,sat_sublattice.py for rule 175 differs from the lines R <= 5 of logs/sat_sublattice_rule175.log)
	@$(H) cd $(W)/rule95 && $(PY) sublattice_any.py 175 5 > out_any175.txt \
	  && diff out_any175.txt <(head -n `wc -l < out_any175.txt` logs/sublattice_any_rule175.log) > /dev/null \
	  || $(call FAIL,sublattice_any.py 175 5 differs from the lines R <= 5 of logs/sublattice_any_rule175.log)
	@$(H) cd $(W)/rule95 && $(PY) sat_sublattice.py 3 all cadical153 250 > out_250.txt \
	  && diff <($(NT) out_250.txt) <($(NT) logs/sat_sublattice_rule250_R3.log) > /dev/null \
	  || $(call FAIL,sat_sublattice.py 3 all cadical153 250 differs from logs/sat_sublattice_rule250_R3.log)

census: setup
	@echo "paper-E census (SAT): centralizer.py --rules 95 --R 1,2 with CaDiCaL and with MiniSat"
	@$(H) cd $(W)/census && $(PY) centralizer.py --rules 95 --R 1,2 --max 5000 | grep -E '^(===|  classification)' > out_cadical.txt \
	  && diff <($(NT) out_cadical.txt) <(grep -E '^(===|  classification)' logs/census_95_R123_cadical153.txt | head -4 | $(NT)) > /dev/null \
	  || $(call FAIL,centralizer.py (cadical153) differs from logs/census_95_R123_cadical153.txt)
	@$(H) cd $(W)/census && $(PY) centralizer.py --rules 95 --R 1,2 --max 5000 --solver minisat22 | grep -E '^(===|  classification)' > out_minisat.txt \
	  && diff <($(NT) out_minisat.txt) <(grep -E '^(===|  classification)' logs/census_95_R123_minisat22.txt | head -4 | $(NT)) > /dev/null \
	  || $(call FAIL,centralizer.py (minisat22) differs from logs/census_95_R123_minisat22.txt)
else
rule95-sat census: setup
	@echo "PySAT (python-sat) not installed: skipping SAT checks of paper-E ($@)"
endif

# ---------------------------------------------------------------- sh
sh: setup
	@echo "paper-E sh: sh 94 0 0 0 0 1 10 12 0"
	@$(H) cd $(W)/sh && $(CC) -O2 -o sh sh.c && ./sh 94 0 0 0 0 1 10 12 0 > out_sh.txt \
	  && test "`grep -o 'objects [0-9]* fail [0-9]*' out_sh.txt | tr '\n' ' '`" = "objects 0 fail 0 objects 1 fail 0 objects 1 fail 0 objects 1 fail 0 objects 2 fail 0 objects 4 fail 0 objects 7 fail 0 objects 12 fail 0 objects 21 fail 0 objects 37 fail 0 " \
	  || $(call FAIL,sh 94 0 0 0 0 1 10 12 0)

# ---------------------------------------------------------------- translation95
translation95: setup
	@echo "paper-E translation95: t95ex.py 9 0, t95ex.py 6 1"
	@$(H) cd $(W)/translation95 && $(CC) -O2 -shared -fPIC -o libgs.so gs.c \
	  && $(PY) t95ex.py 9 0 > out_t95ex.txt && $(PY) t95ex.py 6 1 >> out_t95ex.txt \
	  && sub out_t95ex.txt 'K 9 pairs 3200 failures 0' && sub out_t95ex.txt 'K 6 pairs 200 failures 0'

# ---------------------------------------------------------------- full
full:
	@echo "Long runs of paper E.  They are not executed here; each line gives the"
	@echo "directory, the command, the committed log it reproduces, and its cost."
	@echo
	@echo "code/census (needs python-sat):"
	@echo "  python3 centralizer.py --rules 95 --R 1,2,3 --max 5000 --out census_R3.json   -> logs/census_95_R123_cadical153.txt, ../rule95/census_R3.json  (10 min)"
	@echo "  python3 centralizer.py --rules 95 --R 1,2,3 --max 5000 --solver minisat22     -> logs/census_95_R123_minisat22.txt  (10 min)"
	@echo "  python3 enum_r4.py 95                -> logs/census_95_R4_minisat22.txt   (33 min, 6.6 GB; writes raw_solutions_95_R4.txt, 69 MB)"
	@echo "  python3 enum_r4.py 95 cadical153     -> logs/census_95_R4_cadical153.txt  (47 min, 17.6 GB)"
	@echo
	@echo "code/rule95:"
	@echo "  python3 sat_sublattice.py 7 aut [minisat22]   -> logs/sat_sublattice_aut_R7.log, logs/sat_sublattice_aut_R7_minisat22.log  (25 to 55 min)"
	@echo "  python3 sat_post.py 7 aut                     -> logs/sat_post_R7_aut.log, sat_others_R7_S.npz"
	@echo "  python3 sat_sublattice.py R all|aut minisat22 -> logs/sat_sublattice_all_R1-5_minisat22.log, logs/sat_sublattice_aut_R1-6_minisat22.log"
	@echo "  python3 sat_sublattice.py R all cadical153 175, R = 1..7 (also minisat22) -> logs/sat_sublattice_rule175.log, logs/sat_sublattice_rule175_minisat22.log  (15 s per type at R = 7)"
	@echo "  python3 sublattice_any.py 175 7 [20000 minisat22] -> logs/sublattice_any_rule175.log, logs/sublattice_any_rule175_minisat22.log"
	@echo "  python3 verify_sub.py 5 sat_R5_all.pkl        -> logs/verify_sub_R5.log    (15 min)"
	@echo "  python3 census_cmp.py 4 ../census/raw_solutions_95_R4.txt -> logs/census_cmp_R4.log  (3 min, 1 GB)"
	@echo "  gcc -O2 -fopenmp -o verify7 verify7.c; ./verify7 tx7.bin ty7.bin      -> logs/verify7.log  (about 10 core-minutes)"
	@echo "  ./verify7 tx7.bin ty7.bin 777                                         -> logs/verify7_negative_control.log"
	@echo "  python3 gen36.py                              -> tables36.npz  (2 min)"
	@echo "  python3 enum_golden.py 8                      -> logs/enum_golden_8.log    (25 s)"
	@echo "  python3 enum_golden_linear.py 7               (30 s; expected output in code/rule95/README.md)"
	@echo "  python3 inj_witness.py 10                     (27 s; expected output in code/rule95/README.md)"
	@echo "  python3 sat_sublattice.py 6 aut               -> logs/sat_sublattice_aut_R1-6.log  (15 s)"
	@echo "  python3 block_test.py                         (20 s; expected output in code/rule95/README.md)"
	@echo
	@echo "code/closure:"
	@echo "  ./dc RULE 4, RULE = 91, 94, 155               -> logs/dc_R1-4.log   (about 2 min for 155)"
	@echo "  ./dc 155 5 (or dc_omp)                        -> logs/dc155_R5.log"
	@echo
	@echo "code/rule155:"
	@echo "  python3 mkbase.py 155 6 \"1,1;-1,1\" base6.bin; gcc -O3 -fopenmp -o lc lc.c; ./lc dc 155 6 2 base6.bin -> lc155_R6.log  (90 min with OpenMP on a multi-core machine)"
	@echo
	@echo "code/kink:"
	@echo "  python3 t_R1.py 10                            -> logs/t_R1_10.log   (5 min)"
	@echo
	@echo "code/rule94 (expected output in code/rule94/README.md):"
	@echo "  python3 verify1.py                            (75 s)"
	@echo "  python3 verify2.py 10 12                      (100 s);   python3 verify2.py 8 11 direct"
	@echo "  python3 closure_bf.py 6 13 13 2               (4 s)"
	@echo "  python3 tc.py 12 13; python3 tc.py 16 10; python3 tc.py 8 0 94 1500; python3 tc.py 20 0 94 300  (seconds each)"
	@echo "  python3 te.py 94 11 4,6,8,10,12               (40 s);    python3 te.py 94 12 4,6,8,10,12,13  (2.5 min)"
	@echo
	@echo "code/sh:"
	@echo "  ./sh 94 c0 c1 d0 d1 1 23 24 0 (four tail types), ./sh 91 c0 c1 d0 d1 1 16 20 0 (nine tail types) -> logs/sh_all.log  (minutes to hours per run)"

clean:
	rm -rf .check
