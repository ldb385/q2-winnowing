
import pandas as pd
import os
import qiime2

def summarize( output_dir: str, data: list ) -> None:

    in_featureOrdering, in_AucOrdering, in_PermanovaOrdering = data[0] # Split directory into its parts

    # <><><> Allow user to download data <><><>
    # save out feature ordering data for download
    in_featureOrdering.index.name = 'id'
    in_featureOrdering = qiime2.Metadata( in_featureOrdering.to_frame())
    in_featureOrdering.save(os.path.join( output_dir, 'feature_ordered.tsv'))
    # save out AUC ordering data for download
    in_AucOrdering.index.name = 'id'
    in_AucOrdering = qiime2.Metadata( in_AucOrdering.to_frame())
    in_AucOrdering.save(os.path.join( output_dir, 'auc_ordered.tsv'))
    # save out PERMANOVA ordering data for download
    in_PermanovaOrdering.index.name = 'id'
    in_PermanovaOrdering = qiime2.Metadata( in_PermanovaOrdering.to_frame())
    in_PermanovaOrdering.save(os.path.join( output_dir, 'permanova_ordered.tsv'))

    # index = os.join.path()




    return
