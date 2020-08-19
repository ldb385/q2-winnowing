
import shutil
import pandas as pd

from .plugin_setup import plugin
from ._format import ( WinnowedDirectoryFormat, WinnowedFeatureOrderingFormat,
                       WinnowedAucOrderingFormat, WinnowedPermanovaOrderingFormat )


# <><><> DEFINE STATIC HEADERS FOR SNIFFING <><><>
_EXPECTED_FEATURE_HEADERS_FALSE_ = ["ab_comp", "dataframe1", "metric", "centrality", "total select", "iteration select", "min count",
                    "smooth type", "conditioning", "keep threshold", "correlation", "weighted",
                    "correlation property", "run time", "kappa", "agreement"]
_EXPECTED_FEATURE_HEADERS_TRUE_ = ["ab_comp", "dataframe1", "dataframe2", "metric", "centrality", "iteration select", "total select", "min count",
                    "smooth type", "conditioning", "keep threshold", "correlation", "weighted",
                    "correlation property", "min connected", "run time", "kappa", "agreement"]
_EXPECTED_AUC_HEADERS_ = ["auc", "otu.num"]
_EXPECTED_PERMANOVA_HEADERS_ = ["test", "order", "auc", "SumsOfSqs", "MeanSqs", "F.model", "R2", "Pval",
                                "N.taxa", "F.model.scale"]
    # Placed here since is used in multiple places


# <><><> Define transformers for reading in and writing out format information <><><>
# Define transformer to convert featureOrdering file format into dataframe
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _1( WinnowedFile: WinnowedFeatureOrderingFormat ) -> pd.DataFrame:
    # No Doc String since annotation describes functionality
    with WinnowedFile.open() as wf:
        wf.seek(0)
        wf.readline() # Move past headers
        is_ab_comp = wf.readline().split("\t")[1] # read ab_bool from ab_comp column
        wf.seek(0)  # Reset read cursor to start of data
        if( is_ab_comp == False or is_ab_comp == "False" ):
            expected_col = _EXPECTED_FEATURE_HEADERS_FALSE_ + \
                           list( range(1, abs( len(_EXPECTED_FEATURE_HEADERS_FALSE_) - max([ len(lwf.rstrip("\n").split("\t")) for lwf in wf ]) )) )
            wf.seek( 0 ) # Reset read cursor to start of data
            wf.readline() # skip headers
            orderedFeatures_df = pd.DataFrame()
            for line in wf:
                row = line.rstrip("\n").split("\t")
                if( str(row[0]).isdigit() ):
                    row = row[1:] # index column can be removed
                if( len(row) > len( expected_col ) ):
                    raise ValueError( f"FeatureOrdering: header does not have enough columns for row:\n\t{row}" )
                orderedFeatures_df = orderedFeatures_df.append( [row], ignore_index=True )

            # Make sure dataframe has proper index and column values
            orderedFeatures_df.columns = expected_col # Append changes columns so got to reset
            orderedFeatures_df.reset_index( drop=True, inplace=True )
            return orderedFeatures_df
        else:
            expected_col = _EXPECTED_FEATURE_HEADERS_TRUE_ + \
                           list(range(1, abs(len(_EXPECTED_FEATURE_HEADERS_TRUE_) - max(
                               [len(lwf.rstrip("\n").split("\t")) for lwf in wf]))))
            wf.seek(0)  # Reset read cursor to start of data
            wf.readline()  # skip headers
            orderedFeatures_df = pd.DataFrame()
            for line in wf:
                row = line.rstrip("\n").split("\t")
                if (str(row[0]).isdigit()):
                    row = row[1:]  # index column can be removed
                if (len(row) > len(expected_col)):
                    raise ValueError(f"FeatureOrdering: header does not have enough columns for row:\n\t{row}")
                orderedFeatures_df = orderedFeatures_df.append([row], ignore_index=True)

            # Make sure dataframe has proper index and column values
            orderedFeatures_df.columns = expected_col  # Append changes columns so got to reset
            orderedFeatures_df.reset_index(drop=True, inplace=True)
            return orderedFeatures_df



# Define transformer to convert featureOrdering dataframe to file format
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _2( data: pd.DataFrame ) -> WinnowedFeatureOrderingFormat:
    # No Doc String since annotation describes functionality
    wf = WinnowedFeatureOrderingFormat()
    path = str( wf.path )
    data.to_csv( path, sep="\t" ) # pandas conveniantly converts to tsv

    return wf


# Define transformer to convert auc file format into dataframe
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _3( AucFile: WinnowedAucOrderingFormat ) -> pd.DataFrame:
    # No Doc String since annotation describes functionality
    with AucFile.open() as wf:
        expected_col = _EXPECTED_AUC_HEADERS_
        auc_df = pd.DataFrame()
        wf.readline() # skip headers
        for line in wf:
            row = line.rstrip("\n").split("\t")
            if( str(row[0]).isdigit() ):
                row = row[1:] # index column can be removed
            if( len(row) > len( expected_col ) ):
                raise ValueError( f"AUC: header does not have enough columns for row:\n\t{row}" )
            auc_df = auc_df.append( [row], ignore_index=True )

        # Make sure dataframe has proper index values
        auc_df.columns = expected_col # Append changes columns so got to reset
        auc_df.reset_index( drop=True, inplace=True )
        return auc_df


# Define transformer to convert auc dataframe to file format
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _4( data: pd.DataFrame ) -> WinnowedAucOrderingFormat:
    # No Doc String since annotation describes functionality
    wf = WinnowedAucOrderingFormat()
    path = str( wf.path )
    data.to_csv( path, sep="\t" ) # pandas conveniantly converts to tsv

    return wf


# Define transformer to convert permanova file format into dataframe
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _5( PermanovaFile: WinnowedPermanovaOrderingFormat ) -> pd.DataFrame:
    # No Doc String since annotation describes functionality
    with PermanovaFile.open() as wf:
        expected_col = _EXPECTED_PERMANOVA_HEADERS_
        permanova_df = pd.DataFrame()
        wf.readline() # skip headers
        for line in wf:
            row = line.rstrip("\n").split("\t")
            if( str(row[0]).isdigit() ):
                row = row[1:] # index column can be removed
            if( len(row) > len( expected_col ) ):
                raise ValueError( f"PERMANOVA: header does not have enough columns for row:\n\t{row}" )
            permanova_df = permanova_df.append( [row], ignore_index=True )

        # Make sure dataframe has proper index values
        permanova_df.columns = expected_col # Append changes columns so got to reset
        permanova_df.reset_index( drop=True, inplace=True )
        return permanova_df


# Define transformer to convert permanova dataframe to file format
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _6( data: pd.DataFrame ) -> WinnowedPermanovaOrderingFormat:
    # No Doc String since annotation describes functionality
    wf = WinnowedPermanovaOrderingFormat()
    path = str( wf.path )
    data.to_csv( path, sep="\t" ) # pandas conveniantly converts to tsv

    return wf


# <><><> Define directory transformers for reading in and writing out format information <><><>
# Define transformer to convert permanova file format into dataframe
    # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _7( dirfmt: WinnowedDirectoryFormat ) -> list:
    # No Doc String since annotation describes functionality

        # Ordered feature read of Dir
    artifact_ordered_features = dirfmt.featureOrdering.view( WinnowedFeatureOrderingFormat )
        # PERMANOVA read of Dir
    artifact_permanova = dirfmt.permanova.view( WinnowedPermanovaOrderingFormat )
        # AUC feature read of Dir
    artifact_auc = dirfmt.auc.view( WinnowedAucOrderingFormat )

    dirTuple = [ (
        artifact_ordered_features.view(pd.DataFrame),
        artifact_auc.view(pd.DataFrame),
        artifact_permanova.view(pd.DataFrame)
    ) ]

    return dirTuple


# # Define transformer to convert permanova dataframe to file format
#     # Standard to use a non-meaningful name for plugin transformer
@plugin.register_transformer
def _8( data: list ) -> WinnowedDirectoryFormat:
    # No Doc String since annotation describes functionality
    result = WinnowedDirectoryFormat()
    path = result.path
    of, oa, op = data[0]

    features_fp = str(path / "feature_ordered.tsv")
    permanova_fp = str(path / "permanova_ordered.tsv")
    auc_fp = str(path / "auc_ordered.tsv")

    shutil.copyfile( _2( of ).path, features_fp ) # feature ordering
    shutil.copyfile( _4( oa ).path, auc_fp ) # AUC
    shutil.copyfile( _6( op ).path, permanova_fp ) # PERMANOVA

    return result




