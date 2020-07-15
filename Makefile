
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
    rm -vf !('q2_winnowing/output/$(wildcard *_dump.txt))
    rm -vf !(q2_winnowing/step1_3/output/$(wildcard *_dump.txt)'|'q2_winnowing/step1_3/output/combined_metric_results.csv')
    rm -vf !(q2_winnowing/step4_5/output/$(wildcard *_dump.txt))
    rm -vf !(q2_winnowing/step6/output/$(wildcard *_dump.txt))
    rm -vf !(q2_winnowing/step7_9/output/$(wildcard *_dump.txt))

clean-verbose:
    rm -rf q2_winnowing/output/$(wildcard *_dump.txt)
    rm -rf q2_winnowing/step1_3/output/$(wildcard *_dump.txt)
    rm -rf q2_winnowing/step4_5/output/$(wildcard *_dump.txt)
    rm -rf q2_winnowing/step6/output/$(wildcard *_dump.txt)
    rm -rf q2_winnowing/step7_9/output/$(wildcard *_dump.txt)
