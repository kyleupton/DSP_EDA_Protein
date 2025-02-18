from sklearn.neighbors import kneighbors_graph
from sklearn.cluster import AgglomerativeClustering


def get_AggloModel(key, heatMapData):
    model = False
    if key == 'EuclideanWard_T':
        model = AgglomerativeClustering(n_clusters=None,
                                        metric='euclidean',
                                        memory=None,
                                        connectivity=kneighbors_graph(heatMapData.T, 2),
                                        compute_full_tree=True,
                                        linkage='ward',
                                        distance_threshold=0.1,
                                        compute_distances=True
                                        )
    if key == 'CoseineWard_T':
        model = AgglomerativeClustering(n_clusters=None,
                                        metric='cosine',
                                        memory=None,
                                        connectivity=kneighbors_graph(heatMapData.T, 2),
                                        compute_full_tree=True,
                                        linkage='ward',
                                        distance_threshold=0.1,
                                        compute_distances=True
                                        )
    elif key == 'EuclideanWard_None':
            model = AgglomerativeClustering(n_clusters=None,
                                            metric='euclidean',
                                            memory=None,
                                            connectivity=None,
                                            compute_full_tree=True,
                                            linkage='ward',
                                            distance_threshold=0.1,
                                            compute_distances=True
                                            )

    if key == 'EuclideanWard':
        model = AgglomerativeClustering(n_clusters=None,
                                        metric='euclidean',
                                        memory=None,
                                        connectivity=kneighbors_graph(heatMapData, 2),
                                        compute_full_tree=True,
                                        linkage='ward',
                                        distance_threshold=0.1,
                                        compute_distances=True
                                        )
    return(model)
