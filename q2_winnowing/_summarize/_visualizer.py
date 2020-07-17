
import pandas as pd
import os
import pkg_resources
import qiime2


TEMPLATES = pkg_resources.resource_filename("q2_winnowing", "_summarize")

def summarize( output_dir: str, data: list ) -> None:

    featureOrdering, AucOrdering, PermanovaOrdering = data[0] # Split directory into its parts

    # <><><> Allow user to download data <><><>
    # save out feature ordering data for download
    featureOrdering.index.name = 'id'
    featureOrdering = qiime2.Metadata( featureOrdering.to_frame())
    featureOrdering.save(os.path.join( output_dir, 'feature_ordered.tsv')) # download
    # save out AUC ordering data for download
    AucOrdering.index.name = 'id'
    AucOrdering = qiime2.Metadata( AucOrdering.to_frame())
    AucOrdering.save(os.path.join( output_dir, 'auc_ordered.tsv')) # download
    # save out PERMANOVA ordering data for download
    PermanovaOrdering.index.name = 'id'
    PermanovaOrdering = qiime2.Metadata( PermanovaOrdering.to_frame())
    PermanovaOrdering.save(os.path.join( output_dir, 'permanova_ordered.tsv')) # download

    # Get html files that will be written to
    index_html = os.path.join( TEMPLATES, "assets", "index.html" )
    featureOrdering_html = os.path.join( TEMPLATES, "assets", "featureOrdering.html" )
    auc_html = os.path.join( TEMPLATES, "assets", "auc.html" )
    permanova_html = os.path.join( TEMPLATES, "assets", "permanova.html" )

    # Write tables to specific html files
    featureOrdering.to_html( open( featureOrdering_html, "w" ))
    AucOrdering.to_html( open( auc_html, "w" ))
    PermanovaOrdering.to_html( open( permanova_html, "w" ))


    return
