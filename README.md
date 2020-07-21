# q2-winnowing
<p> 
Qiime2 plugin for inferring the interaction type of microbial communities.</br>
This plugin was developed in order integrate R and Jupyter code into Python as well as Qiime2.</br>
</p>
<b> NOTE THIS IS UNFINISHED </b>
 
# Installation
__NOTE: These instructions assume you are working within a terminal. If not in a terminal please open one to follow along.__
1. Make sure you are running in a qiime2 conda environment or else "import qiime2" will not work.
    * see https://docs.qiime2.org/2020.6/install/ for more information. I recommend using a VM if plugin is being tested or git is being forked.
    * if you have a qiime2 environment created activate with "source activate __env_name__".
        * Ex) if environment is named "qiime2-dev" use "source activate qiime2-dev".
2. clone plugin repo in desired location.
    * run "git clone 'https://github.com/ldb385/q2-winnowing.git'" in desired directory.
3. navigate to inside directory.
    * run "cd q2-winnowing" if git was copied into current working directory ( check with "ls" if unsure if it is in CWD ).
4. "make install"   --__OR__--   "python setup.py install"
    * If error please run python plugin_setup.py.
    * If qiime2 does not import it is likely an issue with PYTHONPATH or conda environment.
    * see qiime2 documentation for more information (https://dev.qiime2.org/latest/tutorials/first-plugin-tutorial/).
5. run "pip install rpy2".
6. run "pip install minepy".
7. run "qiime winnowing".
    * There should be quite a bit of output to console. This is just R configuring libraries and is nothing to be concerned about.
    * If you do not get output it might be that your libraries are already configured. __NOTE: output should only be extensive on first run of command.__
8. Interface should now display in terminal showing __Usage, Options, Commands__, if this does not show please consult previous steps or Qiime2 forums.


# Usage

# Versions
<ul>
 <li> 20.0.0 = Initial Plugin </li>
</ul>
