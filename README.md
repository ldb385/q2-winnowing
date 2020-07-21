# q2-winnowing
<p> 
Qiime2 plugin for inferring the interaction type of microbial communities.</br>
This plugin was developed in order integrate R and Jupyter code into Python as well as Qiime2.</br>
</p>
<b> NOTE THIS IS UNFINISHED </b>
 
# Installation
1. Make sure you are running in a qiime2 conda environment or else "import qiime2" will not work.
    * see https://docs.qiime2.org/2020.6/install/ for more information. I recommend using a VM if plugin is being tested or git is being forked.
    * if you have a qiime2 environment created activate with "source activate __env_name__"
        * if environment is named "qiime2-dev" use "source activate qiime2-dev"
2. pip install rpy2
3. pip install minepy
4. python setup.py install   --OR--   make install
    * If error please run python plugin_setup.py
    * If qiime2 does not import it is likely an issue with PYTHONPATH or conda environment
    * see qiime2 documentation for more information (https://dev.qiime2.org/latest/tutorials/first-plugin-tutorial/)

# Usage

# Versions
<ul>
 <li> 20.0.0 = Initial Plugin </li>
</ul>
