# Top-level targets.
#   make check   quick checks of every paper (minutes on one core)
#   make pdf     compile the five papers (pdflatex, three passes each)
#   make full    print the long runs of every paper and the logs they reproduce
#   make clean   remove scratch directories and LaTeX build products
PYTHON ?= python3
CC ?= gcc
PAPERS = A B C D E

.PHONY: check full pdf clean
check:
	@for p in $(PAPERS); do $(MAKE) --no-print-directory -C paper-$$p check PYTHON="$(PYTHON)" CC="$(CC)" || exit 1; done
	@echo "all quick checks passed"

full:
	@for p in $(PAPERS); do $(MAKE) --no-print-directory -C paper-$$p full PYTHON="$(PYTHON)" CC="$(CC)" || exit 1; done

pdf:
	@for p in $(PAPERS); do \
	  ( cd paper-$$p/tex && for i in 1 2 3; do pdflatex -interaction=nonstopmode -halt-on-error main.tex > /dev/null || { echo "paper-$$p: pdflatex failed, see paper-$$p/tex/main.log"; exit 1; }; done \
	    && echo "paper-$$p/tex/main.pdf" ) || exit 1; done

clean:
	@for p in $(PAPERS); do $(MAKE) --no-print-directory -C paper-$$p clean; done
	@rm -f paper-*/tex/*.aux paper-*/tex/*.log paper-*/tex/*.out paper-*/tex/*.toc
