# Paper C: quick checks and long runs.
#   make check   reruns the quick tier in a scratch copy (.check/) and compares with code/logs/
#   make full    prints the long runs and the logs they reproduce
#   make clean   removes .check/
PYTHON ?= python3
CC ?= gcc
P := $(abspath .)

# $(call run,NAME,COMMAND): run COMMAND in .check/code, compare its output with code/logs/NAME.txt
define run
	@cd .check/code && ( $(2) ) > ../$(1).txt 2> ../$(1).err || { echo "paper-C: $(1) FAILED:"; tail -5 $(P)/.check/$(1).txt $(P)/.check/$(1).err; exit 1; }
	@cmp -s .check/$(1).txt code/logs/$(1).txt || { echo "paper-C: the output of $(1) differs from code/logs/$(1).txt"; diff .check/$(1).txt code/logs/$(1).txt | head -10; exit 1; }
	@echo "  ok  $(1)"
endef
# $(call runsat,NAME,COMMAND,PATTERN): run a SAT program and look for PATTERN in its output (the logs contain timings)
define runsat
	@cd .check/code && ( $(2) ) > ../$(1).txt 2> ../$(1).err || { echo "paper-C: $(1) FAILED:"; tail -5 $(P)/.check/$(1).txt $(P)/.check/$(1).err; exit 1; }
	@grep -q -e $(3) .check/$(1).txt || { echo "paper-C: expected line not found in the output of $(1): "$(3); exit 1; }
	@echo "  ok  $(1)"
endef

.PHONY: check full clean
check:
	@echo "paper-C: quick checks"
	@rm -rf .check && mkdir .check && cp -r code .check/code
	@# --- fixed words (Proposition 5.2(3), Corollary 1.2, Table 1) and walls (Theorems 1.1, 1.3)
	$(call run,fixed_words,$(PYTHON) fixed_words.py)
	$(call run,true_classes,$(PYTHON) true_classes.py)
	$(call run,trains_true,$(PYTHON) trains_true.py)
	$(call run,sym54,$(PYTHON) sym54.py)
	$(call run,cycles54,$(PYTHON) cycles54.py)
	$(call run,walls54,$(PYTHON) walls54.py)
	$(call run,walls54b,$(PYTHON) walls54b.py)
	$(call run,tau_exact,$(PYTHON) tau_exact.py)
	$(call run,verify133,$(PYTHON) verify133.py)
	@# --- Theorem 2.4 and Lemma 2.3
	$(call run,reach_thmA,$(PYTHON) reach_thmA.py)
	$(call run,reach_stair,$(PYTHON) reach_stair.py)
	@# --- hidden walls: rules 59, 27, 127, 179, 223
	$(call run,hw59_walls,$(PYTHON) hw59_walls.py)
	$(call run,hw59_table,$(PYTHON) hw59_table.py)
	$(call run,hw59_overlap,$(PYTHON) hw59_overlap.py)
	$(call run,hw59_classes,$(PYTHON) hw59_classes.py)
	$(call run,hw59_trains,$(PYTHON) hw59_trains.py)
	$(call run,hw59_generators,$(PYTHON) hw59_generators.py)
	$(call run,hw27_walls,$(PYTHON) hw27_walls.py)
	$(call run,hw27_beta,$(PYTHON) hw27_beta.py)
	$(call run,hw27_trains,$(PYTHON) hw27_trains.py)
	$(call run,hwtp_cycles,$(PYTHON) hwtp_cycles.py)
	$(call run,hwtp_walls,$(PYTHON) hwtp_walls.py)
	$(call run,hwtp_trains,$(PYTHON) hwtp_trains.py)
	$(call run,hw_flipall,$(PYTHON) hw_flipall.py 127$(comma)179$(comma)223)
	@# --- rotating walls: rules 219, 126
	$(call run,rot_walls,$(PYTHON) rot_walls.py)
	$(call run,rot_DE,$(PYTHON) rot_DE.py)
	$(call run,rot_hidden,$(PYTHON) rot_hidden.py)
	$(call run,rot_verify,$(PYTHON) rot_verify.py)
	$(call run,rot_trains,$(PYTHON) rot_trains.py)
	@# --- linear and affine rules
	$(call run,markers,$(PYTHON) markers.py)
	$(call run,window_graph,$(PYTHON) window_graph.py)
	$(call run,window_graph_90,$(PYTHON) window_graph.py 90 8)
	$(call run,affine_dim,for f in 170 102 153 150; do $(PYTHON) affine_dim.py $$f 1 2; done)
	$(call run,span_90_R3,$(PYTHON) span90_R3.py)
	$(call run,verify_90_structure,$(PYTHON) verify_90_structure.py)
	$(call runsat,rank_90_R3,$(PYTHON) rank90_packed.py 90 3,"rule 90R R=3: classes 4096$(comma) rank 4031$(comma) kernel dim 65")
	$(call runsat,rank_150_R3,$(PYTHON) rank90_packed.py 150 3,"rule 150R R=3: classes 4096$(comma) rank 4084$(comma) kernel dim 12")
	@# --- SAT censuses (need python-sat)
	@if $(PYTHON) -c "import pysat" 2>/dev/null; then $(MAKE) --no-print-directory check-sat PYTHON="$(PYTHON)"; \
	 else echo "  PySAT (python-sat) not installed: skipping the SAT checks (censuses of 54R at radius 1 and 2, and of 90R at radius 1)"; fi
	@echo "paper-C: all quick checks passed"

comma := ,

.PHONY: check-sat
check-sat:
	$(call runsat,census_54_90_R1,$(PYTHON) centralizer.py --rules 54$(comma)90 --R 1,"rule 54R$(comma) radius R=1: 128 variables$(comma) 17764 clauses.* 13 solutions")
	@grep -q "'named': 5, 'const': 2, 'other': 6, 'invertible': 5" .check/census_54_90_R1.txt || { echo "paper-C: census of 54R at radius 1: unexpected classification"; exit 1; }
	@grep "table=" .check/census_54_90_R1.txt | sort > .check/t1 && grep "table=" code/logs/census_54_90_R1.txt | sort > .check/t2 && cmp -s .check/t1 .check/t2 || { echo "paper-C: the tables of census_54_90_R1 differ from code/logs/census_54_90_R1.txt"; exit 1; }
	$(call runsat,invertible_54_R2,$(PYTHON) invertible_sat.py --rules 54 --R 2 --native,"invertible: 23; of which sigma^a F^b: 13 (expected 13); OTHER invertible: 10")
	@grep "table=" .check/invertible_54_R2.txt | sort > .check/t1 && grep "table=" code/logs/invertible_54_R2.txt | sort > .check/t2 && cmp -s .check/t1 .check/t2 || { echo "paper-C: the tables of invertible_54_R2 differ from code/logs/invertible_54_R2.txt"; exit 1; }

full:
	@echo "paper-C: long runs (not executed by make; run from code/). See CLAIMS.md and code/README.md."
	@echo "  python3 centralizer.py --rules 54 --R 2 --max 20000 --noinj      # logs/census_54_R2_full.txt, 31 min"
	@echo "  python3 invertible_sat2.py --rule 54 --R 3 --max 20001           # logs/balanced_54_R3.txt"
	@echo "  python3 enum54_R3.py run --workers 24; python3 enum54_R3.py classify   # incomplete enumeration, no log"
	@echo "  python3 centralizer.py --rules 59 --R 1,2 --out census_59_R12.json   # logs/census_59_R12.txt, 3 min"
	@echo "  python3 hw59_census.py                                           # logs/hw59_census.txt, 2 min"
	@echo "  python3 centralizer.py --rules 150 --R 1,2 --noinj               # logs/census_150_R12.txt, 30 s"
	@echo "  python3 centralizer.py --rules 90 --R 2 --max 200000 --noinj     # logs/census_90_R2.txt (enumeration, 12 min)"
	@echo "  python3 linear_90.py 1,2                                         # logs/linear_90_R12.txt (first lines)"
	@echo "  python3 affine_dim.py RULE 1 2 3   (RULE = 170, 102, 153, 150)   # logs/affine_dim_R3.txt, 25 s"
	@echo "  python3 hw_search.py 4; python3 hw_search.py 5                   # logs/hw_search_4.txt, logs/hw_search_5.txt"
	@echo "  python3 analyze_90.py 2 all                                      # logs/analyze_90_R2.txt, 1 h 40 min"
	@echo "  python3 tau.py                                                   # logs/tau.txt, 30 s, 4.5 GB"
	@echo "  python3 hw59_period.py                                           # logs/hw59_period.txt, 1 min"
	@echo "  python3 walls_search.py 54 110 30                                # logs/walls_search.txt, 2 min"
	@echo "  python3 hw_search_k6.py                                          # logs/hw_search_k6.txt"
	@echo "  python3 marker_cells.py                                          # logs/marker_cells.txt"
	@echo "  gcc -O2 -o marker_S3 marker_S3.c && ./marker_S3                  # logs/marker_S3.txt"
	@echo "  gcc -O2 -fopenmp -o rot_psi2 rot_psi2.c && ./rot_psi2            # logs/rot_psi2.txt"

clean:
	@rm -rf .check
