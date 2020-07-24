
import pandas as pd
import shutil
import os
import pkg_resources
import qiime2


TEMPLATES = pkg_resources.resource_filename("q2_winnowing", "_summarize")

def summarize( output_dir: str, data: list ) -> None:

    feature_ordering, auc_ordering, permanova_ordering = data[0] # Split directory into its parts

    # <><><> Allow user to download data <><><>
    # save out feature ordering data for download
    feature_ordering_new = os.path.join( output_dir, 'feature_ordered.tsv')
    feature_ordering.to_csv( feature_ordering_new, sep="\t" )
    # save out AUC ordering data for download
    auc_ordering_new = os.path.join( output_dir, 'auc_ordered.tsv')
    auc_ordering.to_csv( auc_ordering_new, sep="\t" )
    # save out PERMANOVA ordering data for download
    permanova_ordering_new = os.path.join( output_dir, 'permanova_ordered.tsv')
    permanova_ordering.to_csv( permanova_ordering_new, sep="\t" )

    # Get html files that will be written to
    index_html = os.path.join( TEMPLATES, "assets", "index.html" )
    feature_ordering_html = os.path.join( TEMPLATES, "assets", "feature_ordering.html" )
    auc_html = os.path.join( TEMPLATES, "assets", "auc.html" )
    permanova_html = os.path.join( TEMPLATES, "assets", "permanova.html" )

    # Write files to output directory
    index_html_new = os.path.join( output_dir, "index.html" )
    shutil.copyfile( index_html, index_html_new )
    feature_ordering_html_new = os.path.join( output_dir, "feature_ordering.html" )
    shutil.copyfile( feature_ordering_html, feature_ordering_html_new )
    auc_html_new = os.path.join( output_dir, "auc.html" )
    shutil.copyfile( auc_html, auc_html_new )
    permanova_html_new = os.path.join( output_dir, "permanova.html" )
    shutil.copyfile( permanova_html, permanova_html_new )

    # Write tables to specific html files
    feature_ordering.to_html( open( feature_ordering_html_new, "w" ))
    auc_ordering.to_html( open( auc_html_new, "w" ))
    permanova_ordering.to_html( open( permanova_html_new, "w" ))

    # copy css stylesheets
    css = os.path.join( TEMPLATES, "assets", "css" )
    css_new = os.path.join( output_dir, "css" )
    shutil.copytree( css, css_new )


    return # Standard to not return on qiime2 visualizer
