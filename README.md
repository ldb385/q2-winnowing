# q2-winnowing
<p> 
Qiime2 plugin for generating a feature selection of data in order to generate a winnowed community.
Diversity and Environmental measures can then be performed on this output to measure the connectivity of the generated winnowed community.</br>
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
    * Other libraries that may need to be installed are:
        * seaborn, matplotlib, sklearn, networkx, scipy, pandas, and numpy
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
               * "qiime tools import --input-path test_otu_table.transpose.biom --type 'FeatureTable\[RelativeFrequency]' --source-format BIOMV210Format --output-path test_otu_table.transpose.qza"
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
        * *qiime winnowing process --i-infile1 __inputFile.qza__ --p-metric graph_centrality --m-sample-types-file __metadata_samples.txt__ --m-sample-types-column Type --p-evaluation kl_divergence --p-min-count 3 --p-conditioning add_one --p-total-select 25 --p-iteration-select {1,4,16,64,128} --p-centrality betweenness --p-keep-threshold 0.5 --p-correlation spearman --p-weighted --p-correlation-prop both --p-detailed --p-verbose --output-dir __./OutFolder__*
        * --__OR__--
        * *qiime winnowing process --i-infile1 __inputFile.qza__ --p-metric graph_centrality --m-sample-types-file __metadata_samples.txt__ --m-sample-types-column Type --p-evaluation kl_divergence --p-min-count 3 --p-conditioning add_one --p-total-select 25 --p-iteration-select {1,4,16,64,128} --p-centrality betweenness --p-keep-threshold 0.5 --p-correlation spearman --p-weighted --p-correlation-prop both --p-detailed --p-verbose --o-result __./OutFile__*

* "qiime winnowing summarize" 
    * is only used to generate a .qzv artifact from the .qza artifact generated from "qiime winnowing processing".
    * To view the .qzv artifact please drag&drop file into https://view.qiime2.org/ .
    * __Command Example)__
        * *qiime winnowing summarize --i-data __./output_file.qza__ --o-visualization __./output_visual.qzv__*

# File Structure Diagram
<!DOCTYPE html>
<html>
<head>
<title>file_structure.html</title>
<meta charset="utf-8"/>
</head>
<body>
<div class="mxgraph" style="max-width:100%;border:1px solid transparent;" data-mxgraph="{&quot;highlight&quot;:&quot;#0000ff&quot;,&quot;nav&quot;:true,&quot;resize&quot;:true,&quot;xml&quot;:&quot;&lt;mxfile host=\&quot;www.diagrameditor.com\&quot; modified=\&quot;2020-08-07T16:52:43.313Z\&quot; agent=\&quot;Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/84.0.4147.105 Safari/537.36\&quot; etag=\&quot;hrjFs0aD30dKfwPm_uIm\&quot; version=\&quot;12.1.3\&quot; type=\&quot;device\&quot; pages=\&quot;1\&quot;&gt;&lt;diagram id=\&quot;2S_ZinG8kKINPOd-_2vE\&quot; name=\&quot;Page-1\&quot;&gt;7Zvdc6M2EMD/Gj/GA4gvP+bycdeZZnptctfmiZFhbasBRISI7fvrK5kPg0SauyuGcaYvMVoJSfx20e4KZYaukt1HhrPNHY0gnllGtJuh65llmbbhix8p2ZcSF6FSsGYkqhodBffkG1RCo5IWJIK805BTGnOSdYUhTVMIeUeGGaPbbrMVjbujZngNmuA+xLEu/ZNEfFNKfcc4yj8BWW/qkU2jqklw3bgS5Bsc0W1LhG5m6IpRysurZHcFsYRXcynvu32ltpkYg5R/zw3rpPjbuv30Fej20va36ObRxxdW2csLjovqgWeWG4v+PuQZTuWs+b5C4T4XcqofthvC4UJUh0J+KZqklCU4PjYQV2v5+2xdbEmaUvFnXfcqpld2XDapyDRjWIwWaQRyxmY91H050vVWGJiQbXgSV9UrEsdXNKbscC+KHPAjW06dM/oErRrfWiLXbcZ7AcZh9ypGs1GOsGqgCXC2F02qG+xKnftucXs0DrPW+KZlGHU7XNnjuun4qDJxUWntBzSINA3e4ScQaGBYuBj8VdgH1w19WK6GgWuaXboLne6iBy46FVxXg/vHzeX13c08iQalu/JDCHvpLn3HFovNSeiaaGK8noY3B15k82w/LN0VuP10I2+xNE5F15uYrq/RfbaC1op8dkuvSrgJEqYibJo97lPBKtx+Ji9jkj51GYqnZvu/RMGYO3XxsV13LZ/baEr7qlQOAZEWuCgcRaSE2Rr4W+5f593i6fTwrGUMYszJS3cafZCrET5TIiZ4VKfRVadvKHrKacFCqO5qRzZKRwvFLFylnxKD1s9B5c1T/wcr6AuizskK/EmtwHHehxXogViQF0mCmcxpztCdeQpPZ+rFVvdntOBZwQeFC6bwZV4f3IXrITyUJ1PhTh3omgsNbs4hMwN0jnGCQhdNHefWCWGLLoec54OyXRqAwO1ja4Bv+P5p2LquMzFbPQbL4mJN0uCccwmFsmd686k561FOmUsMTXiaNcJzpiesRxD1ttmKHiIdfTNOVlzkh11TuRVn+tlO34ire1nWgvl83tqQW6rthKwc8JVNOkGdd1XYVVVKU1D0WolwTNapKIZCSyDkH6QOSYjjy6oiIVEkh+m1maNVDfWaLbqxp+v3GYHdZwTOqYzA/t8IJjaCnmjBdsY0AX3XUaM/fD4JO8Jbt4nSY6vmeJMs/HgOWiZx//LQ6BUljZOEiiy0YwO295NJqNoRGjkLtfQt1bNNkhpPW7OcekPV0jPQICAp4UFwnkGQStiZPFHS09CMZCBWOXgnhKe2YaSnSzLRtwPnDOlqyejUHwSQniS9f9ddL4tvf23wpvTxyod73/pJF6/sM3tqqHBiD4/0LPH2l19vgt++PHz+8nCObl4Bave8xKY76lus52AhTZbCCUWBHIGEAYO8iHk+D/OX94C8Z29/ZOROr1syAxRERZLN+e4sI1iVc882/8ic9dzyMorIYVm1jDpVMOTxoVyOn0az8tAcMFleQwoMc5BSAcwyivxwvd2APCyGQ9GRWKPl6QLjK7AlzWW3M8vFidRH+VdIroFjMURU0sTi2UT/8zPUsON2Ndz3qaGJUtoa9k+mYT0DvBUUBeOUy1N6Rh4ykvGjdm9L7Q5Lf5z4T6Xf93r1wVeT8uHg6/lLfRpvYL6jHMdT+TZr1VuAT3base64BfgjaVasYRmPcyhPZdx3qGlkxnqK+DshCRzQZhCSFQlbPmJ1cAQ58GrdLzL5/IdPcM0aEzKo3cILyQsck2+DLznjfJnT1PW9r8TJ1hxbTwgeIK9ol657WNDjfGhWQdvW6d4LUTyewy8TteN/M6CbfwA=&lt;/diagram&gt;&lt;/mxfile&gt;&quot;,&quot;toolbar&quot;:&quot;pages zoom layers lightbox&quot;,&quot;page&quot;:0}"></div>
<script type="text/javascript" src="https://www.draw.io/js/viewer.min.js"></script>
</body>
</html>


# Versions
<ul>
 <li> 20.0.0 = Initial Plugin </li>
</ul>

# Useful Recourses
* Qiime2 User Documentation: https://docs.qiime2.org/2020.6/
