.PHONY: test test-fast test-scaling verify clean

test: test-fast

test-fast:
	python -m pytest -q -m 'not slow'

test-scaling:
	python -m pytest -q tests/test_v6_scaling.py

verify:
	python -m plasmarom verify --out results_local/verify

clean:
	rm -rf results_local/ .pytest_cache/ build/ dist/
