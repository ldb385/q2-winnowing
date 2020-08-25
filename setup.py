
from setuptools import setup, find_packages

import re
import ast

# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('q2_winnowing/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    __version__ = str(ast.literal_eval(hit))

# setup information of plugin
setup(
    name="winnowing",
    version=__version__,
    packages=find_packages(),
    package_data={"": ["auc.html","auc_ordered.tsv","feature_ordered.html","feature_ordered.tsv",
                       "permanova.html","permanova_ordered.tsv","index.html",
                       "button_formatting.css","dataframe_formatting.css","frame_formatting.css","page_formatting.css"]},
    include_package_data=True,
    author="Kevin Stanley",
    author_email="kgs325@usask.ca",
    description="Used to perform a feature selection on data in order to generate a winnowed community."
                " Diversity and Environmental measures can then be performed on this output to measure the connectivity"
                " of the generated winnowed community.",
    url="https://github.com/ldb385/q2-winnowing",
    entry_points={
        "qiime2.plugins":
            ["q2-winnowing=q2_winnowing.plugin_setup:plugin"]
    },
    zip_safe=False
)
