
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
    rm -rf q2_winnowing/

clean-verbose:
    rm -rf q2_winnowing/output/$(wildcard *_dump.txt)
    rm -rf q2_winnowing/step1_3/$(wildcard *_dump.txt)
    rm -rf q2_winnowing/step4_5/$(wildcard *_dump.txt)
    rm -rf q2_winnowing/step6/$(wildcard *_dump.txt)
    rm -rf q2_winnowing/step7_9/$(wildcard *_dump.txt)
