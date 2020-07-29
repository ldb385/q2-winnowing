
import os
import pandas as pd
import biom

from unittest import TestCase, main as unittest_main
from q2_winnowing.winnow import process

class PluginRunTests( TestCase ):

    package = 'q2_winnowing.tests'
    testing_data = biom.load_table( f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/test_run/test_in_data.biom" )
    testing_metadata = pd.read_csv( f"{os.path.dirname(os.path.realpath(__file__))}/sample_data/test_run/test_in_metadata.txt", sep="\t" ).to_dense()

    def test_plugin_F_gc_kl_ao_cl_sp_bo(self):
        try:
            process( infile1=self.testing_data, ab_comp=False, infile2=None, sample_types=self.testing_metadata,
                     metric="graph_centrality", evaluation="kl_divergence", min_count=3, conditioning="add_one",
                     total_select=25, iteration_select={64,128}, centrality="closeness", keep_threshold=0.5,
                     correlation="spearman", weighted=True, correlation_prop="both", detailed=False, verbose=False)
        except Exception as e:
            self.fail( f"param (False, graph_centrality, kl_divergence, add_one, closeness, spearman, both ) "
                       f"Raised Exception - {e} - unexpectedly.")


    def test_plugin_F_gc_kl_he_ei_sp_ne(self):
        try:
            process( infile1=self.testing_data, ab_comp=False, infile2=None, sample_types=self.testing_metadata,
                     metric="graph_centrality", evaluation="kl_divergence", min_count=3, conditioning="hellinger",
                     total_select=25, iteration_select={64,128}, centrality="eigenvector", keep_threshold=0.5,
                     correlation="spearman", weighted=True, correlation_prop="negative", detailed=False, verbose=False)
        except Exception as e:
            self.fail( f"param (False, graph_centrality, kl_divergence, hellinger, eigenvector, spearman, negative ) "
                       f"Raised Exception - {e} - unexpectedly.")


    def test_plugin_F_pi_kl_ao_bt_pe_po(self):
        try:
            process( infile1=self.testing_data, ab_comp=False, infile2=None, sample_types=self.testing_metadata,
                     metric="pca_importance", evaluation="kl_divergence", min_count=3, conditioning="add_one",
                     total_select=25, iteration_select={64,128}, centrality="betweenness", keep_threshold=0.5,
                     correlation="pearson", weighted=True, correlation_prop="positive", detailed=False, verbose=False)
        except Exception as e:
            self.fail( f"param (False, pca_importance, kl_divergence, add_one, betweenness, pearson, positive ) "
                       f"Raised Exception - {e} - unexpectedly.")


    def test_plugin_F_pi_kl_he_cl_pe_bo(self):
        try:
            process( infile1=self.testing_data, ab_comp=False, infile2=None, sample_types=self.testing_metadata,
                     metric="pca_importance", evaluation="kl_divergence", min_count=3, conditioning="hellinger",
                     total_select=25, iteration_select={64,128}, centrality="degree", keep_threshold=0.5,
                     correlation="pearson", weighted=True, correlation_prop="both", detailed=False, verbose=False)
        except Exception as e:
            self.fail( f"param (False, pca_importance, kl_divergence, hellinger, degree, pearson, both ) "
                       f"Raised Exception - {e} - unexpectedly.")


    def test_plugin_F_lt_kl_he_cl_ke_ne(self):
        try:
            process( infile1=self.testing_data, ab_comp=False, infile2=None, sample_types=self.testing_metadata,
                     metric="log_transform", evaluation="kl_divergence", min_count=3, conditioning="hellinger",
                     total_select=25, iteration_select={64,128}, centrality="closeness", keep_threshold=0.5,
                     correlation="kendall", weighted=True, correlation_prop="negative", detailed=False, verbose=False)
        except Exception as e:
            self.fail( f"param (False, log_transform, kl_divergence, hellinger, closeness, kendall, negative ) "
                       f"Raised Exception - {e} - unexpectedly.")


    def test_plugin_F_ab_kl_he_cl_ke_po(self):
        try:
            process( infile1=self.testing_data, ab_comp=False, infile2=None, sample_types=self.testing_metadata,
                     metric="abundance", evaluation="kl_divergence", min_count=3, conditioning="hellinger",
                     total_select=25, iteration_select={64,128}, centrality="eigenvector", keep_threshold=0.5,
                     correlation="kendall", weighted=True, correlation_prop="positive", detailed=False, verbose=False)
        except Exception as e:
            self.fail( f"param (False, abundance, kl_divergence, hellinger, eigenvector, kendall, positive ) "
                       f"Raised Exception - {e} - unexpectedly.")


    def test_plugin_F_ab_kl_he_cl_mi_bo(self):
        try:
            process( infile1=self.testing_data, ab_comp=False, infile2=None, sample_types=self.testing_metadata,
                     metric="abundance", evaluation="kl_divergence", min_count=3, conditioning="hellinger",
                     total_select=25, iteration_select={64,128}, centrality="betweenness", keep_threshold=0.5,
                     correlation="MIC", weighted=True, correlation_prop="both", detailed=False, verbose=False)
        except Exception as e:
            self.fail( f"param (False, abundance, kl_divergence, hellinger, betweenness, MIC, both ) "
                       f"Raised Exception - {e} - unexpectedly.")



if __name__ == '__main__':
    unittest_main()