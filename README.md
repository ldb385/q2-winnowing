# q2-winnowing
<p> 
Qiime2 plugin for generating a feature selection of data in order to generate a winnowed community.
Diversity and Environmental measures can then be performed on this output to measure the connectivity of the generated winnowed community.</br>
This plugin was developed in order integrate R and Jupyter code into Python as well as Qiime2.</br>
</p>
 
# Installation
__NOTE: These instructions assume you are working within a terminal. If not in a terminal please open one to follow along.__
1. Make sure you are running in a qiime2 conda environment or else "python import qiime2" will not work.
    * See https://docs.qiime2.org/2020.6/install/ for more information. I recommend using a VM if plugin is being tested or git is being forked.
    * If you have a qiime2 environment created activate with "source activate __env_name__".
        * Ex) If environment is named "qiime2-dev" use *source activate qiime2-dev*
1. winnowing packages can be accessed through Conda and Git
    1. Install using Conda. https://anaconda.org/ldb385/q2-winnowing
        1. <pre>conda install -c ldb385 q2-winnowing
        </pre>
    1. Install using Git.
        1. Clone plugin repo in desired location.
            * Run *git clone 'https://github.com/ldb385/q2-winnowing.git'* in desired directory. --__OR__-- Download and unzip zip file in desired directory.
        1. Navigate to inside directory.
            * Run *cd q2-winnowing* if git was copied into current working directory ( check with "ls" if unsure if it is in CWD ).
        1. *make install*   --__OR__--   *python setup.py install*
            * If error please run python plugin_setup.py.
            * If qiime2 does not import it is likely an issue with PYTHONPATH or conda environment.
            * See qiime2 documentation for more information https://dev.qiime2.org/latest/tutorials/first-plugin-tutorial/.
1. <pre>pip install minepy
   pip install rpy2
   pip install seaborn
   </pre>
    * Other libraries used are:
        * matplotlib, sklearn, networkx, scipy, pandas, and numpy
1. Run *qiime winnowing*
    * There should be quite a bit of output to console. This is just R configuring libraries and is nothing to be concerned about.
    * If you do not get output it might be that your libraries are already configured. __NOTE: output should only be extensive on first run of command.__
1. Interface should now display in terminal showing __Usage, Options, Commands__, if this does not show please consult previous steps or Qiime2 forums.

# Usage
__NOTE: These instructions assume you are working within a conda environment in terminal and have completed the Installation portion.__
* Usage should be fairly intuitive however there is some things to note about plugin. Will be easiest to run through potentially confusing parameters for each function. __For examples of files used and produced by the winnowing plugin please see:__<br/>__https://github.com/ldb385/q2-winnowing/tree/master/example_data__
* *qiime winnowing process*
    * --i-infile1 must be frequency table in qza format
       * .qza formatted file can be generated from .csv file. The steps are slightly confusing but it goes .csv-->.txt-->.biom-->.qza .
           1. To convert from .csv it is easiest to save as .txt ( tab delimited ) while in MS Excel.
           2. For next steps a good tutorial is https://cduvallet.github.io/posts/2018/03/qiime2-plugin.
               * In tutorial look for:
               * <pre>biom convert \
                    -i test_otu_table.transpose.txt \
                    -o test_otu_table.transpose.biom \
                    --table-type="OTU table" \
                    --to-hdf5
                 </pre>
               * <pre>qiime tools import \
                    --input-path test_otu_table.transpose.biom \
                    --type 'FeatureTable[Frequency]' \
                    --input-format BIOMV210Format \
                    --output-path test_otu_table.transpose.qza
                  </pre>
       * The same follows when using --i-infile2.
    * --m-sample-types-file
       * This is the path to the .txt metadata ( tab delimited ). Format is:
       * <pre>#SampleID    Type
         sample1      natural
         sample2      natural
         sample3      natural
           ...          ...
         sample15     invaded
         sample16     invaded
         sample17     natural
           ...          ...
         </pre>
    * --m-sample-types-column "Type"
    * __Command Example)__
        * <pre>qiime winnowing process \
             --i-infile1 <b>inputFile.qza</b> \
             --p-metric graph_centrality \
             --m-sample-types-file <b>metadata_samples.txt</b> \
             --m-sample-types-column Type \
             --p-evaluation kl_divergence \
             --p-min-count 3 \
             --p-conditioning add_one \
             --p-total-select 25 \
             --p-iteration-select {1,4,16,64,128} \
             --p-centrality betweenness \
             --p-keep-threshold 0.5 \
             --p-correlation spearman \
             --p-weighted \
             --p-correlation-prop both \
             --p-detailed \
             --output-dir <b>./OutFolder</b> \
             --verbose
          </pre>
        * --__OR__--
        * <pre>qiime winnowing process \
             --i-infile1 <b>inputFile.qza</b> \
             --p-metric graph_centrality \
             --m-sample-types-file <b>metadata_samples.txt</b> \
             --m-sample-types-column Type \
             --p-evaluation kl_divergence \
             --p-min-count 3 \
             --p-conditioning add_one \
             --p-total-select 25 \
             --p-iteration-select {1,4,16,64,128} \
             --p-centrality betweenness \
             --p-keep-threshold 0.5 \
             --p-correlation spearman \
             --p-weighted \
             --p-correlation-prop both \
             --p-detailed \
             --o-result <b>./OutFile</b> \
             --verbose
          </pre>

* *qiime winnowing summarize* 
    * is only used to generate a .qzv artifact from the .qza artifact generated from "qiime winnowing processing".
    * To view the .qzv artifact please drag&drop file into https://view.qiime2.org/ .
    * __Command Example)__
        * <pre>qiime winnowing summarize \
             --i-data <b>./output_file.qza</b> \
             --o-visualization <b>./output_visual.qzv</b> \
             --verbose
          </pre>

# File Structure Diagram
![q2-winnowing-plugin-structure](https://user-images.githubusercontent.com/55117132/91479374-9144da80-e85e-11ea-9f62-f704b210efb2.png)

# Versions
<ul>
 <li> 20.0.0 = Initial Plugin </li>
</ul>

# Useful Resources
* Qiime2 User Documentation: https://docs.qiime2.org/2020.6/
* Biom Format Documentation: https://biom-format.org/documentation/biom_conversion.html
* Anaconda Documentation: https://docs.anaconda.com/anaconda-cloud/user-guide/tasks/work-with-packages/

