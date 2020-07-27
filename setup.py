
from setuptools import setup, find_packages

# NOTE: This is temporary placement and is subject to change
# Nothing specified should be interpreted as factual or unchangeable

import re
import ast

# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('q2_winnowing/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    __version__ = str(ast.literal_eval(hit))

# setup information of plugin
# TODO: Change to valid information
setup(
    name="winnowing",
    version=__version__,
    packages=find_packages(),
    author="",
    author_email="",
    description="Infer the interaction type of microbial communities through statistical analysis. " +
                "This will allow for a better understanding of taxa interaction at a micro scale.",
    url="https://github.com/ldb385/q2-winnowing",
    entry_points={
        "qiime2.plugins":
            ["q2-winnowing=q2_winnowing.plugin_setup:plugin"]
    },
    zip_safe=False
)