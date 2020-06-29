# Loading packages here
import sys
import argparse
import numpy as np
# import skbio as sb  # Use this for ANOSIM
import pandas as pd  # Use this for working with dataframes
import os
import networkx
import scipy.stats as stats
# import scipy.spatial.distance as distance # conda install scipy
# from skbio.diversity import beta_diversity # for bray curtis conditioning
from sklearn.decomposition import PCA  # Use for PCA
import matplotlib.pyplot as plt  # Use this for plotting
# from skbio.stats.distance import anosim
import math
import csv
import minepy  # pip install minepy
import time
import seaborn as sb

global min_count
global window_size
global outdir
global disjoint
global connectedness


# Verbose global will help prevent not necessary graphs from always being generated
global verbose


def remove_min_count(df, min_count):
    """
    Function remove_min_count: This function removes data that is all zeros in a column
        best used once merging has taken place to get rid of all features that are zero in both conditions
    :param df: @type pandas dataframe: The data to remove counts below min_count
    :return: @type pandas dataframe: The resulting dataframe after removal
    """
    return (df.loc[:, (df > min_count).any(axis=0)])


def add_one_smoothing(df):
    """
    Function add_one_smoothing: Add one accounts for the possibility of unseen events (0s) occurring in the future.
    :param df: @type pandas dataframe: The data to smooth
    :return: @type pandas dataframe: The smoothed data
    """
    temp = df.copy() + 1
    temp = temp / temp.sum()
    return temp


def bray_curtis(df):
    """
    Function bray_curtis: Performs a bray-curtis dissimilarity
    :param df: @type pandas dataframe: The data to do the dissimilarity on
    :return: @type pandas dataframe: The bray-curtis'ed data
    Computes the Bray-Curtis distance between two 1-D arrays.
    """
    temp = df.copy()
    ids = df.columns.values.tolist()
    bc = beta_diversity("braycurtis", temp.transpose(), ids)
    bc_df = pd.DataFrame(data=bc.data, index=bc.ids, columns=bc.ids)
    return bc_df


def hellinger(df):
    """
    Function hellinger: The hellinger transformation deals with the double zero problem in ecology.
            The hellinger transformation is the square root of the result of the current row
            divided by the sum of all rows. This is done for each element in the dataframe.
    :param df: @type pandas dataframe: The data to do the transformation on.
    :return: @type pandas dataframe: A dataframe of the hellinger transformed data.
    """
    temp = df.copy()
    hellinger_data = np.sqrt(temp.div(temp.sum(axis=1), axis=0))  # Hellinger transformation
    return hellinger_data


def condition(df, cond_type):
    """
    Function condition: A preprocessing step to condition the data based on the type specidied
    :param data: @type pandas dataframe - The data to be conditioned
    :param cond_type: @type string - The type of conditioning to run. Valid values: add_one, hellinger
    :return: @type pandas dataframe - The conditioned dataframe.
    """
    # print('conditioning type', cond_type)
    temp = df.copy()
    temp = temp.loc[(temp != 0).any(axis=1)]
    if cond_type == 'add_one':
        conditioned_data = add_one_smoothing(temp)
    elif cond_type == 'hellinger':
        conditioned_data = hellinger(temp)
    elif cond_type == 'bray_curtis':
        conditioned_data = bray_curtis(temp)
    else:
        conditioned_data = temp
    return conditioned_data


def smooth(df, type):
    """
    Function smooth: Smoothing function to filter out noise
    :param df: @type pandas dataframe: The data to be smoothed
    :param type: @type string: The type of smoothing to do (currently sliding_window is the only option)
    :return: @type pandas dataframe: A dataframe of the smoothed data
    """
    temp = df.copy()
    if type == 'sliding_window':
        result = temp.rolling(window_size, min_periods=1, center=True).mean()
    else:
        result = temp
    return result


def pca_abundance(df, num_pca_components=4, cond_type='hellinger'):
    """
    Function pca_abundance: running PCA iteratively on the data, removing the highest/lowest abundance features
        (based on sorting of abundances array). The gradient is used to find the greatest change in inertia when
         features are removed. Then we select all the features up to that point. These are the most important features
    :param data: @type pandas dataframe: The data to run pca on
    :param num_pca_components: @type integer: The number of components to use for pca
    :param cond_type: @type string: The conditioning type to use for pca (eg. hellinger, add_one)
    :param smoothing_type: @type string: The type of smoothing to do on the dataframe of total eigenvalues found by PCA
    :return: important_features - a list of the most important features as found by running the PCA
    """

    pca = PCA(n_components=num_pca_components)  # Run a PCA with n components
    data = df.copy()  # copy the data into a dataframe to manipulate
    eigen_df = pd.DataFrame()  # create a dataframe to hold the eigenvalues from the pca

    abundances_arr = pd.unique(data.sum(axis=0))  # get the set of unique values
    abundances_arr = np.sort(abundances_arr)[::-1]  # sort highest to lowest

    # now we want to loop through all the abundances and find those with the most variance
    # once we find those we'll have to match them back up to the features with those abundances
    # so that we can send back the list of sorted features
    for i in range(len(abundances_arr)):
        if len(abundances_arr) - i == num_pca_components:
            break
        conditioned_df = condition(data, cond_type)  # condition data
        result = pca.fit(conditioned_df)  # Run the PCA on the conditioned data here
        variance_arr = result.explained_variance_ratio_  # Output the variance associated with each eigenvector
        drop_list = list(
            data.columns[data.sum(axis=0) == abundances_arr[i]])  # Find all features with the current abundance
        components = result.components_
        variance_df = pd.DataFrame(variance_arr, columns=[
            str(abundances_arr[i])]).transpose()  # Convert the eigenvalues to a data frame of OTU rows x N components
        # variance_df = pd.DataFrame(variance_arr).transpose()  # Convert the eigenvalues to a data frame of OTU rows x N components
        eigen_df = eigen_df.append(variance_df)  # Append to the eigenvalue df
        data.drop(drop_list, inplace=True, axis=1)  # Drop all the features with the current abundance
        # You can only iterate over the number of features minus the number of components.
        if len(abundances_arr) - i == num_pca_components:
            break

    eigen_df['Total'] = eigen_df.sum(axis=1)  # sum up the eigenvalues to get the total variance of all components
    # print('eigen df', eigen_df)
    total_eigen = eigen_df.copy().iloc[:, [-1]]
    total_eigen.sort_values(by='Total', ascending=0,
                            inplace=True)  # order the values in descending order, since we want to remove the highest eigenvalues

    # loop through each row and get the feature name and the abundance.
    # Match the feature name to the eigenvector variance from the total_eigen dataframe.
    # Then we'll return a dataframe with the feature as an index and the variance as the value
    ordered_features = pd.DataFrame(columns=['variance'])
    for index, row in df.sum(axis=0).iteritems():
        if str(row) in total_eigen.index:
            # print('variance = ',total_eigen['Total'].loc[str(row)])
            ordered_features.loc[index] = total_eigen['Total'].loc[str(row)]
        ordered_features.sort_values(by='variance', ascending=0, inplace=True)
    return ordered_features


def pca_importance(df, num_pca_components=4, cond_type='hellinger'):
    """
    Function pca_importance: running PCA and selecting the most important features as found
        by the eigenvectors.
    :param data: @type pandas dataframe: The data to run pca on
    :param num_pca_components: @type integer: The number of components to use for pca
    :param select_n: @type integer: The number of important features to return
    :param cond_type: @type string: The conditioning type to use for pca (eg. hellinger, add_one)
    :param smoothing_type: @type string: The type of smoothing to do on the dataframe of total eigenvalues found by PCA
    :return: ordered_features - a numpy array of the features ordered from most important to least, as found by running the PCA
    """

    pca = PCA(n_components=num_pca_components)  # Run a PCA with n components
    data = df.copy()  # copy the data into a dataframe to manipulate
    conditioned_df = condition(data, cond_type)  # condition data

    result = pca.fit(conditioned_df)  # Run the PCA on the conditioned data here
    components = result.components_  # the eigenvectors/components
    eigenvectors = pd.DataFrame(components, columns=[data.columns])

    abs_eigens = np.absolute(eigenvectors)  # element wise absolute value to get rid of negatives
    variance = pd.DataFrame({'metric': abs_eigens.sum(
        axis=0)})  # sum up the components to get the amount of variance across all components for each feature
    variance['pca1'], variance['pca2'] = eigenvectors.iloc[0], eigenvectors.iloc[1]

    # order the features from largest to smallest to return as our sorted dataframe
    ordered_features = variance.sort_values(by='metric', ascending=0)
    # print('ordered_features columns', ordered_features.index.values)
    return ordered_features


def pca_legendre(df, num_pca_components, cond_type='hellinger'):
    return 0


def abundance(df):
    data = df.copy()
    summed_vals = data.sum(axis=0)
    sorted_abundances = summed_vals.sort_values(ascending=0)
    return sorted_abundances


def find_correlation(df, corr_type='spearman'):
    df_r = 0
    data = df.copy()
    if corr_type == 'MIC':
        # the pstats output for the mic and tic are condensed 1D arrays
        # we need to turn the output into a 2D upper triangular matrix and then mirror it to get the
        # full correlation matrix
        micp, ticp = minepy.pstats(data.T, alpha=0.6, c=15, est="mic_approx")
        num_features = data.shape[1]
        tri = np.zeros((num_features, num_features))
        tri[np.triu_indices(num_features, 1)] = micp
        full_corr = tri + tri.T

        df_r = pd.DataFrame(full_corr)
        df_r.columns = data.columns.values
        df_r.index = data.columns.values
    else:
        df_r = data.corr(corr_type)

    if isinstance(df_r, pd.DataFrame):
        df_r.fillna(0, inplace=True)  # ugly hack to make the NAs go away, should work for sampling but not advisable
        df_r = df_r[(df_r != 0).any(axis=1)]
        df_r = df_r.loc[:, (df_r != 0).any(axis=0)]
    return df_r


# this function returns the sorted centrality for a given centrality
# given a dataframe organized as an adjacency matrix, build a graph and compute the centrality
# return sorted centrality and the graph in networkx format
def graph_centrality(df, cent_type='betweenness', keep_thresh=0.5, cond_type='add_one', corr_type='spearman',
                     weighted=False, corr_dir='none', min_connected=0):
    """
    :param df: @type pandas DataFrame
    :param cent_type: @type string - valid values: betweenness, degree, closeness, eigenvector
    :param keep_thresh: @type float - default 0.5
    :param cond_type: @type: string - valid values: add_one, hellinger, bray_curtis
    :param corr_type: @type: string - valid values: spearman, kendall, pearson, MIC
    :param weighted: @type: boolean - True if you want to produce a graph with weighted edges, False otherwise
    :param corr_dir: @type: string - valid values: none, positive, negative
    :return:
    """

    print('In Graph Centrality Function')

    data = df.copy()
    conditioned_df = condition(data, cond_type)  # condition data
    w_corr_df = find_correlation(conditioned_df, corr_type)
    if corr_dir == 'positive':
        w_corr_df_b = 1 - w_corr_df.copy()  # only keep strong positive correlations (small positive numbers)
    elif corr_dir == 'negative':
        w_corr_df_b = 1 + w_corr_df.copy()  # only keep strong negative correlations (small negative numbers)
    else:
        w_corr_df_b = 1 - abs(w_corr_df.copy())  # keep both strong positive and negative correlations
    w_corr_df_b[
        (w_corr_df_b >= 1 - keep_thresh)] = 1  # set anything greater than the threshold value to 1 so we can remove it.
    labels = list(w_corr_df_b.index)
    temp = abs(w_corr_df_b.copy())
    temp.insert(0, 'var1', labels)

    if weighted == True:
        attr = 'weight'
    else:
        attr = 'edge'

    df_b = pd.melt(temp, 'var1', var_name='var2', value_name=attr)

    df_b = df_b.loc[((df_b[attr] <= 1 - keep_thresh) & (df_b[attr] >= 0.0)),
           :]  # take only those edge pairs that are at or above the threshold
    df_g = networkx.from_pandas_edgelist(df_b, 'var1', 'var2', attr)  # takes a list of valid edges

    # create a graphml file
    # graph_filename = "graph-{}.graphml".format(process_id)
    # networkx.write_graphml(df_g, os.path.join(outdir,graph_filename))
    # draw and display the graph with labels
    # networkx.draw(df_g, with_labels=True)
    # networkx.draw(df_g)
    # pylab.show()

    # store the adjacency matrix in a pandas dataframe and then export it to a csv
    # if there isn't an adjacency matrix csv with this id, create a csv file

    # am = networkx.to_pandas_adjacency(df_g)
    # adjacency_filename = "adj_matrix-{}.csv".format(process_id)
    # am.to_csv(os.path.join(outdir,adjacency_filename))

    # check to see if the graph is disjoint. If it is, we just want to return an empty dataframe
    # and stop the winnowing process

    num_subgraphs = networkx.number_connected_components(df_g)
    total_nodes = networkx.number_of_nodes(df_g)
    largest_subgraph = max( (df_g.subgraph(c) for c in networkx.connected_components( df_g ) ), key=len )
    # largest_subgraph = max( networkx.connected_component_subgraphs(df_g), key=len ) #Deprecated
    nodes_sg = networkx.number_of_nodes(largest_subgraph)

    print('total', total_nodes, ' largest', nodes_sg)
    percent_connected = nodes_sg / total_nodes * 100
    print('percent connected', percent_connected)
    global connectedness
    connectedness.append((nodes_sg, total_nodes, float(str(round(percent_connected, 2)))))
    print('connectedness', connectedness)
    if percent_connected < float(min_connected):
        print('graph under min connectedness... returning')
        disjoint = True
        return pd.DataFrame()
    else:
        if cent_type == 'betweenness':
            centrality = networkx.betweenness_centrality(largest_subgraph)
        elif cent_type == 'degree':
            centrality = networkx.degree_centrality(largest_subgraph)
        elif cent_type == 'closeness':
            centrality = networkx.closeness_centrality(largest_subgraph)
        elif cent_type == 'eigenvector':
            try:
                centrality = networkx.eigenvector_centrality(largest_subgraph)
            except:
                print('eigenvector failed to converge... returning')
                disjoint = True
                return pd.DataFrame()
        else:
            # print('error, unknown centrality')
            return -1

        centrality_df = pd.DataFrame.from_dict(centrality, orient='index')
        centrality_df.columns = ['metric']

        if not centrality_df.empty:
            centrality_df = centrality_df[centrality_df.iloc[:, 0] > 0]

        if not centrality_df.empty:
            centrality_df.sort_values('metric', axis=0, ascending=False, inplace=True)

        '''fig = plt.figure()
        plt.hist(centrality_df, bins=20)
        plt.xlabel('Centrality')
        plt.ylabel('Frequency')
        plt.title('Graph Centrality Distribution')
        plt.tight_layout()
        #fig.savefig(os.path.join(outdir,'test.jpg'))
        plt.show()'''

    return centrality_df


def selection(func, s_total, s_per_iter, df, *args):
    """
    Function selection: does N loops of the metric before going to the evaluation step
    :param type: are we selecting features to remove or to retain
    :param N: the number of features for selection
    :return: not sure yet
    """

    data = df.copy()
    feature_list = []
    selected_df = pd.DataFrame()
    # the metric (func) returns a dataframe of the most important features
    # when we get back this dataframe, we need to select a specified number of features (select_per_iter)
    # until we select the total number of features desired (select_total)
    selected = 0

    if s_total == 'all':
        select_total = len(data.columns)
    else:
        try:
            select_total = int(s_total)
        except:
            print('not a valid select total value')
            sys.exit(1)

    if s_per_iter == 'all':
        select_per_iter = len(data.columns)
    else:
        try:
            select_per_iter = int(s_per_iter)
        except:
            print('not a valid select per iteration value')

    # make sure they aren't trying to select more features per iter than total features
    select_per_iter = min(select_per_iter, select_total)
    select = select_per_iter
    for i in range(0, math.ceil(select_total / select_per_iter)):
        # call the metric with the current data
        sorted_df = func(data, *args)

        if not sorted_df.empty:
            if ((i + 1) * select_per_iter > select_total):
                select = select_total % selected
            # take the top n features returned by the metric
            top_features = sorted_df.iloc[:select].index.values

            selected_df = selected_df.append(sorted_df.iloc[:select])
            selected += select

            # if this is the last time we're selecting features,
            # look at the metric value of last feature selected
            # and see if there were others in the list with the same
            # value. Include those too.
            if selected == select_total:
                last_row = selected_df.tail(1)
                last_feature = last_row.index.values[0]
                last_value = last_row.iloc[0]['metric']
                same_vals = sorted_df[(sorted_df['metric'] == last_value) & (~sorted_df.index.isin(selected_df.index))]
                selected_df = selected_df.append(same_vals)

            # add to the list of features selected
            feature_list.extend(top_features)
            #remove the top features from the data frame
            data.drop(top_features.tolist(), axis=1, inplace=True)
        else:
            return selected_df

    results = selected_df
    return results


def reduction(func, select_total, remove_per_iter, df, *args):
    return 0


def evaluation(func, *args):
    result = func(*args)
    return result


def pca_inertia_eval(df, num_pca_components, cond_type, smoothing_type):
    """
    Function pca_abundance: running PCA iteratively on the data, removing the highest/lowest abundance features
        (based on sorting of abundances array). The gradient is used to find the greatest change in inertia when
         features are removed. Then we select all the features up to that point. These are the most important features
    :param data: @type pandas dataframe: The data to run pca on
    :param num_pca_components: @type integer: The number of components to use for pca
    :param cond_type: @type string: The conditioning type to use for pca (eg. hellinger, add_one, bray_curtis)
    :param smoothing_type: @type string: The type of smoothing to do on the dataframe of total eigenvalues found by PCA
    :return: important_features - a list of the most important features as found by running the PCA
    """

    pca = PCA(n_components=num_pca_components)  # Run a PCA with n components
    data = df.copy()  # copy the data into a dataframe to manipulate
    eigen_df = pd.DataFrame()  # create a dataframe to hold the eigenvalues from the pca

    abundances_arr = pd.unique(data.sum(axis=0))  # get the set of unique values

    #   ------------ QUESTION ----------
    # if it's sorted from lowest to highest, it can't run the PCA on the four biggest abundances
    # so they aren't in the eigen_df dataframe, meaning they aren't added to the important OTU list
    # Should we always sort highest to lowest, or should I be creating the important_OTUS list differently?
    abundances_arr = np.sort(abundances_arr)[::-1]  # sort highest to lowest

    for i in range(len(abundances_arr)):
        conditioned_df = condition(data, cond_type)  # condition data
        result = pca.fit(conditioned_df)  # Run the PCA on the conditioned data here
        variance_arr = result.explained_variance_ratio_  # Output the variance associated with each eigenvector
        drop_list = list(data.columns[data.sum(axis=0) == abundances_arr[i]])  # Find the top features
        variance_df = pd.DataFrame(variance_arr, columns=[
            str(abundances_arr[i])]).transpose()  # Convert the eigenvalues to a data frame of OTU rows x N components
        eigen_df = eigen_df.append(variance_df)  # Append to the eigenvalue df
        data.drop(drop_list, inplace=True, axis=1)  # Drop all the features with the current abundance

        # You can only iterate over the number of features minus the number of components.
        if len(abundances_arr) - i == num_pca_components:
            break

    eigen_df['Total'] = eigen_df.sum(axis=1)  # sum up the eigenvalues to get the total variance of all components
    total_eigen = eigen_df.copy().iloc[:, [-1]]
    # sorted list to return
    total_eigen.sort_values(by='Total', ascending=0,
                            inplace=True)  # order the values in descending order, since we want to remove the highest eigenvalues

    # evaluation
    smoothed = smooth(total_eigen, smoothing_type)  # then smooth the values using a sliding window

    # use the gradient function to find the greatest slope between two abundances
    # this is our max inertia and we will return features above that point
    gradient = abs(np.gradient(smoothed.values.flatten()))
    smoothed['Gradient'] = gradient
    maxGradient = smoothed['Gradient'].argmax()
    topAbundances = smoothed.loc[:maxGradient].index.values

    # get all the features that have the abundance totals selected
    important_features = []
    for i in range(len(topAbundances)):
        important_features.extend(data.columns[data.sum(axis=0) == int(topAbundances[i])].values)

    # show the graph of the smoothed summed eigenvalues from running the PCA at each step.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    smoothed['Total'].plot(kind='bar', title='Total Variation')
    ax.set_xticklabels([])
    # plt.show()

    return important_features


def kl_divergence(A, B, features, select_per_iter, cond_type):
    """
    :param A: @type: pandas dataframe - data file 1
    :param B: @type: pandas dataframe
    :param selected_features: @type - list
    :param select_per_iter: @type - int
    :return: @type list - the list of kl divergence values
    """
    selected_features = list(features)
    # condition(data, cond_type)
    data1 = condition(A.sum(axis=0).transpose(),
                      cond_type)  # KL divergence requires that histograms are the same size so sum to remove differences in number of samples
    data2 = condition(B.sum(axis=0).transpose(), cond_type)
    num_features = len(selected_features)
    diverge_vals = []
    feature_count = 0
    select = select_per_iter
    for i in range(0, math.ceil(num_features / select_per_iter)):
        if ((i + 1) * select_per_iter > num_features):
            select = num_features % feature_count

        # need to drop the first 'select' features from the list
        selections = selected_features[:select]
        del selected_features[:select]
        data1.drop(selections, inplace=True)
        data2.drop(selections, inplace=True)

        if len(data1.shape) > 1:  # if there is more than one dimension, flatten
            tempA = data1.values.flatten()
        else:
            tempA = data1.values

        if len(data2.shape) > 1:  # if there is more than one dimension, flatten
            tempB = data2.values.flatten()
        else:
            tempB = data2.values

        feature_count += select
        kl_diverge = stats.entropy(tempA, tempB, 2.0)  # find the KL-Divergence base 2
        if kl_diverge > 1e50:
            kl_diverge = 1e50
        diverge_vals.append(kl_diverge)  # determine if the remaining histograms are more alike after elimination

    return diverge_vals


def create_pca_plot(features):
    data = features.copy()
    plt.scatter(data['pca1'], data['pca2'])
    plt.xlabel('Principle Component 1')
    plt.ylabel('Principle Component 2')
    plt.title('Scatter Plot of Principle Components 1 and 2')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "pca_scatter.png"))
    # plt.show()


def plot_graph_centrality(features, cond_type, corr_type, corr_dir, keep_thresh, weighted):
    data = features.copy()
    # create a graph of the edges/nodes using the same centrality type as used to select the features
    # this is the top25 file stuff

    conditioned_df = condition(data, cond_type)  # condition data
    w_corr_df = find_correlation(conditioned_df, corr_type)

    if corr_dir == 'positive':
        w_corr_df_b = 1 - w_corr_df.copy()  # only keep strong positive correlations (small positive numbers)
    elif corr_dir == 'negative':
        w_corr_df_b = 1 + w_corr_df.copy()  # only keep strong negative correlations (small negative numbers)
    else:
        w_corr_df_b = 1 - abs(w_corr_df.copy())  # keep both strong positive and negative correlations
    w_corr_df_b[
        (w_corr_df_b >= 1 - keep_thresh)] = 1  # set anything greater than the threshold value to 1 so we can remove it.

    labels = list(w_corr_df_b.index)
    temp = abs(w_corr_df_b.copy())
    temp.insert(0, 'var1', labels)

    if weighted == True:
        attr = 'weight'
    else:
        attr = 'edge'

    df_b = pd.melt(temp, 'var1', var_name='var2', value_name=attr)
    df_b = df_b.loc[((df_b[attr] <= 1 - keep_thresh) & (df_b[attr] > 0.0)),
           :]  # take only those edge pairs that made the cut
    df_g = networkx.from_pandas_edgelist(df_b, 'var1', 'var2', attr)  # takes a list of valid edges
    networkx.write_graphml(df_g, 'graph_network.graphml')
    networkx.draw(df_g, node_color='dodgerblue', edge_color='dimgrey', with_labels=True)

    plt.savefig(os.path.join(outdir, "graph_network.png"))
    create_ecological_network(df_b.copy(), data.copy())
    # plt.show()


def create_ecological_network(df, data):
    feature_names = data.columns
    metric_matrix = pd.DataFrame(columns=feature_names, index=feature_names, dtype='float64')
    metric_matrix.fillna(value=0.0, inplace=True)

    for i in df.index:
        row = df.loc[i].tolist()
        metric_matrix.loc[row[0]][row[1]] = row[2]
        # print("Rows:0 "+str(row[0])+" Row:1 "+str(row[1])+" Row:2 "+str(row[2])+" Set Value: "+str(metric_matrix[row[0]][row[1]]))

    metric_matrix.to_csv(os.path.join(outdir, "Metric Network.csv"))

    if (verbose):
        plt.figure()
        plt.title("Feature Metric Heatmap")
        plt.tight_layout()
        sb.heatmap(metric_matrix, annot=True, cmap=['Grey', 'Blue'], cbar=False)
        plt.savefig(fname=os.path.join(outdir, "Metrics Network.png"), format='png', dpi=600, bbox_inches='tight',
                    papertype='ledger')


def plot_feature_metric(features):
    # x-axis is the OTU (feature) in ranked order
    # y-axis is the metric value
    data = features.copy()
    data['metric'].plot(style='.-')
    plt.xticks(np.arange(0, len(data['metric'])), data.index.values, rotation=45, ha='center')
    plt.xlabel('Feature Label')
    plt.ylabel('Metric Value')
    plt.title('Metric Value Per Feature')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "metric_value.png"))
    # plt.show()


def name_imp_features(feature_df, filename):
    """
    Function name_imp_features: takes the dataframe of important features in the format
                                [feature_id,metric] and a file to match the names of
                                the features to the actual names from a given file.
                                This is for naming species.
    :param feature_df: the dataframe of features we want to name
    :param filename: the filename (including absolute path if necessary) of the
                     file to use to map the names of the features. This file must be in
                     the format: feature name, [size], kingdom, phylum, class, order,
                     family, genus, species. Each of those names must be wrapped in
                     a formatter (e.g. k__Bacteria(100)
    :return: the named feature df
    """
    # read in the naming file
    namefile_path = filename
    name_file = pd.read_csv(namefile_path, index_col='OTU')
    # clean data file

    if 'Size' in name_file.columns:
        del name_file['Size']

    name_file.replace(to_replace='[a-z]__', value='', inplace=True, regex=True)
    name_file.replace(to_replace='\([0-9]*\)', value='', inplace=True, regex=True)

    # make the dictionary
    dicty = name_file.T.to_dict('list')
    name_file['name'] = [[i for i in dicty[k] if i][-1] for k in dicty]
    name_df = pd.DataFrame(index=name_file.index)
    name_df['name'] = name_file['name']

    # join the feature id to the actual name
    merged = pd.concat([feature_df, name_df['name']], axis=1, join_axes=[feature_df.index])
    name_dict = merged['name'].to_dict()
    renamed_df = merged.rename(index=name_dict)
    return renamed_df['metric']


def name_feature_list(feature_df, filename):
    """
    Function name_feature_list: takes the dataframe of the winnowed data in its original
                                abundance format and a file to match the names of
                                the features to the actual names from a given file.
                                This is for naming species.
    :param feature_df: the dataframe of features we want to name
    :param filename: the filename (including absolute path if necessary) of the
                     file to use to map the names of the features. This file must be in
                     the format: feature name, [size], kingdom, phylum, class, order,
                     family, genus, species. Each of those names must be wrapped in
                     a formatter (e.g. k__Bacteria(100)
    :return: the named feature df
    """
    # read in the naming file
    namefile_path = filename
    name_file = pd.read_csv(namefile_path, index_col='OTU')
    # clean data file

    if 'Size' in name_file.columns:
        del name_file['Size']

    name_file.replace(to_replace='[a-z]__', value='', inplace=True, regex=True)
    name_file.replace(to_replace='\([0-9]*\)', value='', inplace=True, regex=True)

    # make the dictionary
    dicty = name_file.T.to_dict('list')
    name_file['name'] = [[i for i in dicty[k] if i][-1] for k in dicty]
    name_df = pd.DataFrame(index=name_file.index)
    name_df['name'] = name_file['name']

    # join the feature id to the actual name
    name_dict = name_df['name'].to_dict()
    renamed_df = feature_df.rename(columns=name_dict)
    return renamed_df


def log_transfrom_balance(df, cond_type='add_one'):
    print("In the Log Transfom Balance Function")

    data = condition(df.copy(), cond_type)

    mertic_result = pd.DataFrame(columns=['OTU', 'metric'])

    for col in data.columns:
        k_pos = int(data[col].count() / 2)
        k_neg = int(data[col].count() / 2)
        # data[col]+=1
        sum_k_pos_log = 0
        sum_k_neg_log = 0

        for j in range(0, k_pos):
            sum_k_pos_log += math.log(data[col][j])
            sum_k_neg_log += math.log(data[col][j + k_pos])

        balance = 1 / k_pos * sum_k_pos_log - 1 / k_neg * sum_k_neg_log
        mertic_result = mertic_result.append({'OTU': col, 'metric': balance}, ignore_index=True)

    mertic_result.set_index('OTU', inplace=True)
    return mertic_result


def log_transfrom(df, cond_type='add_one'):
    print("In the log Transfom Function")

    return np.log(condition(df.copy(), cond_type))


def main(ab_comp, dataframe1, dataframe2, metric_name, c_type, min_count,
         total_select, iteration_select, pca_components, smooth_type,
         window_size, centrality_type, keep_threshold, correlation,
         weighted, corr_prop, evaluation_type, plot_metric,
         create_graph, plot_pca, naming_file, proc_id, min_connected,
         detailed=False, verbose_p=False):

    t_start = time.perf_counter()

    global verbose # Using global since this will not change and passing verbose to all methods will be confusing
    verbose = verbose_p

    # read the file(s) into a pandas dataframe and condition it
    dataframe1.reset_index( drop=True, inplace=True )
    dataframe1.fillna(0, inplace=True)

    global disjoint
    disjoint = False

    global process_id
    process_id = proc_id

    global outdir
    outdir = f"{os.path.dirname(os.path.realpath(__file__))}/output/{metric_name}_{correlation}_{str(keep_threshold)}_{centrality_type}"
    resultsOutdir = f"{os.path.dirname(os.path.realpath(__file__))}/output"
    # allows for cleaner execution and use of relative paths
    os.makedirs(outdir, exist_ok=True)

    global connectedness
    connectedness = []

    # set up the metric variables and the file names
    if ab_comp:
        metric_header = ["ab_comp", "dataframe1", "dataframe2", "metric", "centrality", "total select", "iteration select",
                         "min count", "smooth type", "conditioning", "keep threshold", "correlation",
                         "weighted", "correlation property", "min connected", "run time" ]

        metric_params = [ab_comp, dataframe1.name, dataframe2.name, metric_name, centrality_type, total_select,
                         iteration_select, min_count, smooth_type, c_type, keep_threshold, correlation,
                         weighted, corr_prop, min_connected]
        eval_params = [ab_comp, dataframe1.name, dataframe2.name, evaluation, c_type, total_select, iteration_select]

        dataframe2.fillna(0, inplace=True)

        data_file = pd.concat([dataframe2, dataframe1])

        if( detailed ):
            metric_filename = f"{dataframe1.name}-{dataframe2.name}-results-{process_id}.csv"
            abundance_filename = f"{dataframe1.name}-{dataframe2.name}-abundances-{process_id}.csv"

    else:
        metric_header = ["ab_comp", "dataframe1", "metric", "centrality", "total select", "iteration select",
                         "min count", "smooth type", "conditioning", "keep threshold", "correlation",
                         "weighted", "correlation property", "run time" ]

        metric_params = [ab_comp, dataframe1.name, metric_name, centrality_type, total_select,
                         iteration_select, min_count, smooth_type, c_type, keep_threshold, correlation, weighted,
                         corr_prop]
        eval_params = [ab_comp, dataframe1.name, evaluation, c_type, total_select, iteration_select]

        data_file = dataframe1

        if( detailed ):
            metric_filename = f"{dataframe1.name}-{process_id}.csv"
            abundance_filename = f"{dataframe1.name}-abundances-{process_id}.csv"

    if min_count != -1:
        data = remove_min_count(data_file, min_count)
    else:
        data = data_file

    # log_transfrom(data)

    # run the metric selection step to return the important features
    important_features = pd.DataFrame()

    if metric_name == 'graph_centrality':
        metric = graph_centrality
        important_features = selection(metric, total_select, iteration_select, data, centrality_type, keep_threshold,
                                       c_type, correlation, weighted, corr_prop, min_connected)
    elif metric_name == 'pca_importance':
        metric = pca_importance
        important_features = selection(metric, total_select, iteration_select, data, pca_components, c_type)
    elif metric_name == 'abundance':
        metric = abundance
        important_features = reduction(metric, total_select, iteration_select, data)
    elif metric_name == 'log_transform':
        metric = graph_centrality
        important_features = selection(metric, total_select, iteration_select, log_transfrom(data, c_type),
                                       centrality_type, keep_threshold, c_type, correlation, weighted, corr_prop,
                                       min_connected)

    # print("Printing Log Transformed Important Features")
    # print(important_features)

    t_end = time.perf_counter()
    runtime = t_end - t_start
    metric_params.append(runtime)

    # add the abundance totals to the resulting dataframe and create a list of the important feature names
    important_features['abundances'] = data[important_features.index.values].sum(axis=0) #add abundances to the df
    important_feature_list = list(important_features.index.values)

    #create a dataframe with the abundances of the features determined as 'important'
    feature_abundances = data[important_feature_list]

    if important_features.empty:
        important_features.loc[1] = 'Warning: No features returned.'
        if disjoint:
            important_features.loc[2] = 'Graph is disjoint'
    elif naming_file:
        important_features = name_imp_features(important_features, naming_file)
        feature_abundances = name_feature_list(feature_abundances, naming_file)

    if plot_metric:
        plot_feature_metric(important_features)

    if plot_pca and (metric == pca_importance or metric == pca_abundance):
        create_pca_plot(important_features)
    if create_graph and metric == graph_centrality:
        plot_graph_centrality(feature_abundances, c_type, correlation, corr_prop, keep_threshold, weighted)

    #
    # create csv files for all the outputs.
    #

    # add the metric information to the important features and append to the metric results dataframe
    feature_row = metric_params + important_feature_list

    if( detailed ):
        with open( os.path.join( outdir, 'metric_results.csv'), 'a') as f:
            writer = csv.writer(f)
            writer.writerow(feature_row)
        with open( os.path.join( resultsOutdir, 'combined_metric_results.csv'), 'a' ) as f:
            writer = csv.writer(f)
            writer.writerow(feature_row)

    metric_header.extend( range( 1, abs( len(feature_row) - len(metric_header) ) +1 ))
    feature_df = pd.DataFrame( [feature_row], columns=metric_header) # format list into data frame before output

    # get the abundances for each of the important features and write those to a new file
    # print('final important features', important_features)
    if( detailed ):
        important_features.to_csv(os.path.join(outdir, metric_filename))
        feature_abundances.to_csv(os.path.join(outdir, abundance_filename), index=False)

    if( ab_comp ):
        param_dict = {'ab_comp': ab_comp, 'dataframe1': dataframe1.name, 'dataframe2': dataframe2.name,
                      'metric_name': metric_name, 'c_type': c_type, 'min_count': min_count,
                      'total_select': total_select, 'iteration_select': iteration_select,
                      'pca_components': pca_components, 'smooth_type': smooth_type,
                      'window_size': window_size, 'centrality_type': centrality_type,
                      'keep_threshold': keep_threshold, 'correlation': correlation,
                      'weighted': weighted, 'corr_prop': corr_prop, 'evaluation_type': evaluation_type,
                      'plot_metric': plot_metric, 'create_graph': create_graph, 'plot_pca': plot_pca,
                      'naming_file': naming_file, 'proc_id': proc_id, 'runtime': runtime, 'min_connected': min_connected,
                      'connectedness': connectedness}
    else:
        param_dict = {'ab_comp': ab_comp, 'dataframe1': dataframe1.name, 'dataframe2': "",
                      'metric_name': metric_name, 'c_type': c_type, 'min_count': min_count,
                      'total_select': total_select, 'iteration_select': iteration_select,
                      'pca_components': pca_components, 'smooth_type': smooth_type,
                      'window_size': window_size, 'centrality_type': centrality_type,
                      'keep_threshold': keep_threshold, 'correlation': correlation,
                      'weighted': weighted, 'corr_prop': corr_prop, 'evaluation_type': evaluation_type,
                      'plot_metric': plot_metric, 'create_graph': create_graph, 'plot_pca': plot_pca,
                      'naming_file': naming_file, 'proc_id': proc_id, 'runtime': runtime, 'min_connected': min_connected,
                      'connectedness': connectedness}

    if( detailed ):
        parameter_df = pd.DataFrame(list(param_dict.items()),
                                    columns=['Parameter', 'Values'])
        param_filename = f'parameter_list-{process_id}.csv'
        parameter_df.to_csv(os.path.join(outdir, param_filename))

    return ( feature_df, important_features, feature_abundances )

# graph information loss, kl outlier divergence.
# % of total information loss since removal has occurred. or look at inflection point and how far away your N otus are from that.

# front end takes in a file parses each line to set the selection parameters.
# for each selection parameter, put the output files into a folder.
# Then output the folder at the end.
