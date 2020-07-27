# q2-winnowing
<p> 
Qiime2 plugin for inferring the interaction type of microbial communities.</br>
This plugin was developed in order integrate R and Jupyter code into Python as well as Qiime2.</br>
</p>
<b> NOTE THIS IS UNFINISHED </b>
 
# Installation
__NOTE: These instructions assume you are working within a terminal. If not in a terminal please open one to follow along.__
1. Make sure you are running in a qiime2 conda environment or else "python import qiime2" will not work.
    * See https://docs.qiime2.org/2020.6/install/ for more information. I recommend using a VM if plugin is being tested or git is being forked.
    * If you have a qiime2 environment created activate with "source activate __env_name__".
        * Ex) If environment is named "qiime2-dev" use "source activate qiime2-dev".
1. Clone plugin repo in desired location.
    * Run "git clone 'https://github.com/ldb385/q2-winnowing.git'" in desired directory. --__OR__-- Download and unzip zip file in desired directory.
1. Navigate to inside directory.
    * Run "cd q2-winnowing" if git was copied into current working directory ( check with "ls" if unsure if it is in CWD ).
1. "make install"   --__OR__--   "python setup.py install"
    * If error please run python plugin_setup.py.
    * If qiime2 does not import it is likely an issue with PYTHONPATH or conda environment.
    * See qiime2 documentation for more information https://dev.qiime2.org/latest/tutorials/first-plugin-tutorial/.
1. Run "pip install minepy".
1. Run "qiime winnowing".
    * There should be quite a bit of output to console. This is just R configuring libraries and is nothing to be concerned about.
    * If you do not get output it might be that your libraries are already configured. __NOTE: output should only be extensive on first run of command.__
1. Interface should now display in terminal showing __Usage, Options, Commands__, if this does not show please consult previous steps or Qiime2 forums.

# Usage
__NOTE: These instructions assume you are working within a conda environment in terminal and have completed the Installation portion.__
* Usage should be fairly intuitive however there is some things to note about plugin. Will be easiest to run through potentially confusing parameters for each function. 
* "qiime winnowing process"
    * --i-infile1 must be frequency table in qza format
       * .qza formatted file can be generated from .csv file. The steps are slightly confusing but it goes .csv-->.txt-->.biom-->.qza .
           1. To convert from .csv it is easiest to save as .txt ( tab delimited ) while in MS Excel.
           2. For next steps a good tutorial is https://cduvallet.github.io/posts/2018/03/qiime2-plugin.
               * In tutorial look for:
               * "biom convert -i test_otu_table.transpose.txt -o test_otu_table.transpose.biom --table-type="OTU table" --to-hdf5"
               * "qiime tools import --input-path test_otu_table.transpose.biom --type 'FeatureTable[RelativeFrequency]' --source-format BIOMV210Format --output-path test_otu_table.transpose.qza"
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
        * *qiime winnowing process --i-infile1 __inputFile.qza__ --p-metric graph_centrality --m-sample-types-file metadata_samples.txt --m-sample-types-column Type --p-evaluation kl_divergence --p-min-count 3 --p-conditioning add_one --p-total-select 25 --p-iteration-select {1,4,16,64,128} --p-centrality betweenness --p-keep-threshold 0.5 --p-correlation spearman --p-weighted --p-correlation-prop both --p-detailed --p-verbose --output-dir __./OutFolder__*
        * --__OR__--
        * *qiime winnowing process --i-infile1 __inputFile.qza__ --p-metric graph_centrality --m-sample-types-file metadata_samples.txt --m-sample-types-column Type --p-evaluation kl_divergence --p-min-count 3 --p-conditioning add_one --p-total-select 25 --p-iteration-select {1,4,16,64,128} --p-centrality betweenness --p-keep-threshold 0.5 --p-correlation spearman --p-weighted --p-correlation-prop both --p-detailed --p-verbose --o-result __./OutFile__*

* "qiime winnowing summarize" 
    * is only used to generate a .qzv artifact from the .qza artifact generated from "qiime winnowing processing".
    * To view the .qzv artifact please drag&drop file into https://view.qiime2.org/ .
    * __Command Example)__
        * *qiime winnowing summarize --i-data __./output_file.qza__ --o-visualization __./output_visual.qzv__*

# Versions
<ul>
 <li> 20.0.0 = Initial Plugin </li>
</ul>

# Useful Recourses
* Qiime2 User Documentation: https://docs.qiime2.org/2020.6/
