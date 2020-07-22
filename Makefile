
PYTHON ?= python


test:
	py.test

test-cov:
	py.test --cov=q2_winnowing

install:
	$(PYTHON) setup.py install

dev:
	pip install -e .

clean-detailed:
	ls ./q2_winnowing/output/ | grep -v '$(wildcard ./q2_winnowing/output/*_dump.txt)' | xargs rm -f
	ls ./q2_winnowing/step1_3/output/ | grep -v '$(wildcard ./q2_winnowing/step1_3/output/*_dump.txt)' | grep -v 'combined_metric_results.csv' | xargs rm -f
	ls ./q2_winnowing/step4_5/output/ | grep -v '$(wildcard ./q2_winnowing/step4_5/output/*_dump.txt)' | xargs rm -f
	ls ./q2_winnowing/step6/output/ | grep -v '$(wildcard ./q2_winnowing/step6/output/*_dump.txt)' | xargs rm -f
	ls ./q2_winnowing/step7_9/output/ | grep -v '$(wildcard ./q2_winnowing/step7_9/output/*_dump.txt)' | xargs rm -f

clean-verbose:
	rm -f '$(wildcard q2_winnowing/output/*_dump.txt)'
	rm -f '$(wildcard q2_winnowing/step1_3/output/*_dump.txt)'
	rm -f '$(wildcard q2_winnowing/step4_5/output/*_dump.txt)'
	rm -f '$(wildcard q2_winnowing/step6/output/*_dump.txt)'
	rm -f '$(wildcard q2_winnowing/step7_9/output/*_dump.txt)'
