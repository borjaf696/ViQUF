//
// Borja :)
//
#include "clusterContainer.h"

vector<vector<int>> ClusterContainer::getClusters()
{
    return _clusters;
}

void ClusterContainer::classifyClusters()
{
    vector<set<string>> environmentsPerCluster, readsPerCluster;
    _hybridClusters = vector<bool>(_clusters.size(), false);
    size_t currCluster = 0;
    for (auto cluster: _clusters)
    {
        set<string> e_c, e_r;
        for (auto c: cluster) {
            e_c.emplace(_f_environments[_rs_ds(c)]);
            e_r.emplace(_reads[_rs_ds(c)]);
        }
        environmentsPerCluster.push_back(e_c);
        readsPerCluster.push_back(e_r);
        _hybridClusters[currCluster++] = (e_c.size() > 1);
    }
    /*
     * CheckResults
     */
    for (auto set_:environmentsPerCluster)
    {
        cout <<"NewCluster: "<<endl;
        for (auto e:set_)
            cout << e << ",";
        cout <<endl;
    }
    for (auto set_:readsPerCluster)
    {
        cout <<"NewCluster: "<<endl;
        for (auto e:set_)
            cout << e << ",";
        cout <<endl;
    }
}
