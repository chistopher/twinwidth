
#pragma once

#include "utils.hpp"

using Graph = vector<vector<int>>; // adj matrix; 0: no edge 1: edge 2: red edge -1: deleted vertex
struct EdgeList {
    int n;
    vector<pair<int,int>> edges;
    vector<int> idx; // index of original node
};

EdgeList readInput();
vector<EdgeList> components(const EdgeList& g);
Graph toMat(const EdgeList& g);
Graph subgraph(const Graph& g, const vector<int>& nodes);
Graph aliveSubgraph(const Graph& g);
void draw(const Graph& g, map<int,string> color = {}, map<int,string> label= {});
void draw_snapshot(const Graph& g, int time, int v, int p);