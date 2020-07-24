
import pandas as pd
import shutil
import os
import pkg_resources
import qiime2


TEMPLATES = pkg_resources.resource_filename("q2_winnowing", "_summarize")

def summarize( output_dir: str, data: list ) -> None:

    featureOrdering, AucOrdering, PermanovaOrdering = data[0] # Split directory into its parts

    # <><><> Allow user to download data <><><>
    # save out feature ordering data for download
    featureOrdering_new = os.path.join( output_dir, 'feature_ordered.tsv')
    featureOrdering.to_csv( featureOrdering_new, sep="\t" )
    # save out AUC ordering data for download
    AucOrdering_new = os.path.join( output_dir, 'auc_ordered.tsv')
    AucOrdering.to_csv( AucOrdering_new, sep="\t" )
    # save out PERMANOVA ordering data for download
    PermanovaOrdering_new = os.path.join( output_dir, 'permanova_ordered.tsv')
    PermanovaOrdering.to_csv( PermanovaOrdering_new, sep="\t" )

    # Get html files that will be written to
    index_html = os.path.join( TEMPLATES, "assets", "index.html" )
    featureOrdering_html = os.path.join( TEMPLATES, "assets", "feature_ordering.html" )
    auc_html = os.path.join( TEMPLATES, "assets", "auc.html" )
    permanova_html = os.path.join( TEMPLATES, "assets", "permanova.html" )

    # Write files to output directory
    index_html_new = os.path.join( output_dir, "index.html" )
    shutil.copyfile( index_html, index_html_new )
    featureOrdering_html_new = os.path.join( output_dir, "feature_ordering.html" )
    shutil.copyfile( featureOrdering_html, featureOrdering_html_new )
    auc_html_new = os.path.join( output_dir, "auc.html" )
    shutil.copyfile( auc_html, auc_html_new )
    permanova_html_new = os.path.join( output_dir, "permanova.html" )
    shutil.copyfile( permanova_html, permanova_html_new )

    # Write tables to specific html files
    featureOrdering.to_html( open( featureOrdering_html_new, "w" ))
    AucOrdering.to_html( open( auc_html_new, "w" ))
    PermanovaOrdering.to_html( open( permanova_html_new, "w" ))

    # copy css stylesheets
    css = os.path.join( TEMPLATES, "assets", "css" )
    css_new = os.path.join( output_dir, "css" )
    shutil.copytree( css, css_new )


    return # Standard to not return on qiime2 visualizer
