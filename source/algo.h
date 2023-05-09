
#pragma once
#include "graph.h"

struct Solution {
    Graph mat;
    struct TreeRepresentation {
        vector<int> order;
        vector<int> parent;
        vector<int> get_partition() const; // for each node, current living representative
    };
    TreeRepresentation merges;
    int width = 1e9;
    vector<int> rdeg;

    struct MergeData {
        int width_before;
        vector<int> rdeg, adj_from, adj_to;
    };
    vector<MergeData> merge_data = {};

    void merge(int from, int to);
    void pop_merge();
};

struct Algo {
    bool verbose = false;
    map<vector<int>,int> cache;
    Solution best;
    int m_lower_bound = 0;
    int search_space = 0;
    map<vector<int>, vector<pair<int,int>>> cluster;

    Solution solve(const Graph& g, int lb = 0, int ub = 1e9);

private:
    void iterate_trees(Solution& partial); // runtime: fac(n) * fac(n-1) // 2**(n-1) * poly(n)
    void add_cluster(const Solution& partial, int p);
    vector<pair<int,int>> check_clusters(const Solution& partial);
};
