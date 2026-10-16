# Recipe lines must start with a TAB character, not spaces.
PY ?= python
SMOKE_DIR ?= /tmp/codon-smoke
NOTEBOOKS = notebooks/01_wrangling.ipynb notebooks/02_eda.ipynb notebooks/03_preprocessing.ipynb \
            notebooks/04_modeling.ipynb notebooks/05_results.ipynb
# each notebook runs top to bottom in its own fresh kernel, in the order listed
EXECUTE = $(PY) -m jupyter nbconvert --to notebook --execute \
          --ExecutePreprocessor.kernel_name=python3 --ExecutePreprocessor.timeout=-1 \
          --ExecutePreprocessor.record_timing=False

.PHONY: help test smoke run check

help:
	@echo "make test    unit and data tests (about 10 s)"
	@echo "make smoke   all notebooks, tiny design (CODON_SMOKE=1), into $(SMOKE_DIR) (about 1-2 min)"
	@echo "make run     all notebooks in place, full design (about 40 min on 2 cores)"
	@echo "make check   size, errors, stderr and execution order of the committed notebooks"

test:
	$(PY) -m pytest -q tests

smoke:
	mkdir -p $(SMOKE_DIR)
	CODON_SMOKE=1 $(EXECUTE) --output-dir $(SMOKE_DIR) $(NOTEBOOKS)
	$(PY) scripts/check_notebooks.py $(SMOKE_DIR)/*.ipynb

run:
	$(EXECUTE) --inplace $(NOTEBOOKS)
	$(PY) scripts/check_notebooks.py $(NOTEBOOKS)

check:
	$(PY) scripts/check_notebooks.py $(NOTEBOOKS)
