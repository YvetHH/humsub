import numpy as np

import numpy as np

from scipy.spatial import distance


def get_prediction_strength_(closest_centroid, test_labels):
    '''
    Function for calculating the prediction strength of clustering
    
    Parameters
    ----------
    closest_centroid: array
        N x M matrix for 
    test_labels : array
        Labels predicted for the test set
        
    Returns
    -------
    prediction_strength : float
        Calculated prediction strength
    '''
    
    
    n_test = len(test_labels)
    unique_labels = np.unique(test_labels)
    assert closest_centroid.shape[0]== n_test
    
    #assert set(np.unique(closest_centroid)) == set(unique_labels), f"Test labels {unique_labels} != {}
    




    # populate the co-membership matrix
    closest_centroid= np.matrix(closest_centroid)
    Is_same_cluster= ((closest_centroid - closest_centroid.T) == 0 ).astype(int)


    # calculate the prediction strengths for each cluster

    strength_clusters=[]

    for k in unique_labels:

        ids_of_cluster = (test_labels == k) 
        cluster_size= sum(ids_of_cluster)

        # if only one element in cluster inf prediction strength
        if cluster_size <= 1:
            strength_clusters.append(float('inf')) 

        else:
            
            # sum of co-localized , substract self count
            n_colocalized =Is_same_cluster[ids_of_cluster][:,ids_of_cluster].sum()- cluster_size
            # fraction of co-localized elments
            ps= n_colocalized / (cluster_size* (cluster_size-1))
            strength_clusters.append(ps)
            



    # return minimum
    return min(strength_clusters)


def calculate_prediction_strength(test_data,train_centers,test_labels,metric='euclidean'):
    


    DM= distance.cdist(test_data,train_centers,metric=metric)

    # calculate closest train centroid for each test sample
    closest_centroid = np.argmin( DM ,axis=1)
    
    return get_prediction_strength_(closest_centroid,test_labels)


def calculate_prediction_strength_medioid(train_centroids_ids, test_sample_ids, test_labels,distance_matrix):
    

    DM= np.array(distance_matrix)[test_sample_ids][:,train_centroids_ids]


    # calculate closest train centroid for each test sample
    closest_centroid = np.argmin( DM ,axis=1)


    return get_prediction_strength_(closest_centroid,test_labels)




# get median of cluster
import warnings 

def get_medioid_indices(distance_matrix, labels):

    D= np.array(distance_matrix)
    cluster_lables= np.unique(labels)
    
    median_indexes = np.zeros(cluster_lables.shape,dtype=int)
    
    # Update the medoids for each cluster
    for i,k in enumerate(cluster_lables):
        # Extract the distance matrix between the data points
        # inside the cluster k
        cluster_k_idxs = np.where(labels == k)[0]

        if len(cluster_k_idxs) == 0:
            warnings.warn(f"Cluster {k} is empty! ")
            continue

        in_cluster_distances = D[
            cluster_k_idxs, cluster_k_idxs[:, np.newaxis]
        ]

        # Calculate all costs from each point to all others in the cluster
        in_cluster_all_costs = np.sum(in_cluster_distances, axis=1)

        min_cost_idx = cluster_k_idxs[np.argmin(in_cluster_all_costs)]

        median_indexes[i] = min_cost_idx

    return median_indexes


def calculate_clutering( indexes,distance_matrix,k  ):

    distance_matrix_sub = distance_matrix.iloc[indexes,indexes]
    linkage_sub= average( sp.distance.squareform(distance_matrix_sub))
    
    labels = hc.fcluster(linkage_sub,k,criterion='maxclust')
    return labels






def calculate_prediction_strength_labels(labels_from_test,labels_from_train,kmax,verbose=False):

    strength_clusters= []

    used_test_labels=set()

    unique_train_labels, counts_train_labels= np.unique(labels_from_train, return_counts=True)
    
    
    for train_cluster in unique_train_labels[np.argsort(counts_train_labels)]:

        corresponding_test_labels = labels_from_test[labels_from_train==train_cluster]
        
        unique_test_labels, counts_test_labels= np.unique(corresponding_test_labels, return_counts=True)
        
        is_unused_label= ~ np.isin(unique_test_labels, used_test_labels)
        
        

        if any(is_unused_label):
            
            N_matching = counts_test_labels[is_unused_label].max()
            best_matching_label = unique_test_labels[ counts_test_labels==N_matching  ][0]
            
            used_test_labels.add(best_matching_label)
            
            if verbose:
                print(f"Label {train_cluster} in train corresponds to {best_matching_label} in test")
            
            if (N_matching!= counts_test_labels.max()):
                print('suboptiomal choice')
            
            
            # fraction of co-localized elments
            ps= N_matching / corresponding_test_labels.shape[0]
            strength_clusters.append(ps)          
            
        else:
            strength_clusters.append(float('inf')) 
            print(f"Didn't found label for {train_cluster}")






    if verbose:
        print(strength_clusters)
    return min(strength_clusters)

