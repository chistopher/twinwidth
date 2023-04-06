
#include <bits/stdc++.h>

#define rep(a, b)   for(int a = 0; a < (b); ++a)
#define all(a)      (a).begin(),(a).end()
//#define endl        '\n'

using namespace std;
using ll = long long;
using Graph = vector<vector<int>>; // adj matrix; 0: no edge 1: edge 2: red edge -1: deleted vertex

Graph generateER(int n, double p) {
    uniform_real_distribution dist;
    mt19937 gen;
    vector mat(n, vector(n, 0));
    rep(i,n) rep(j,n)
            if(i<j && dist(gen)<p)
                mat[i][j] = mat[j][i] = 1;
    return mat;
}

struct TreeRepresentation {
    vector<int> order;
    vector<int> parent;
};

struct Solution {
    Graph mat;
    TreeRepresentation merges;
    int width = 1e9;
    void merge(int from, int to) {
        rep(i,ssize(mat)) {
            if(i==from || i==to) continue;
            auto same = (mat[i][from] == mat[i][to]);
            auto new_val = same ? mat[i][to] : 2;
            mat[i][to] = mat[to][i] = new_val;

        }
        rep(i,ssize(mat)) mat[i][from] = mat[from][i] = -1;
        rep(i,ssize(mat)) width = max(width, (int)count(all(mat[i]),2));
        merges.order.push_back(from);
        merges.parent[from] = to;
    }
};

optional<Solution> partitionDP(const Graph& g, int upper_bound = 1e9) {
    using Partition = vector<int>; // rep for each merged comp is the highest vertex inside; partition is rep for each vertex
    auto n = ssize(g);

    map<Partition, Solution> cache;
    queue<Partition> q;
    Partition start(n);
    iota(all(start),0);
    cache[start] = {g,{.order = {},.parent = vector(n,-1)},0};
    q.push(start);
    while(size(q)) { // basically a BFS
        Partition curPar = q.front(); q.pop();
        auto remVerts = curPar;
        sort(all(remVerts));
        remVerts.erase(unique(all(remVerts)),end(remVerts)); // now contains all reps of partition
        auto curSol = cache[curPar];
        if(size(remVerts)==1)
            return curSol; // all nodes merged

        for(auto v : remVerts) {
            for(auto u : remVerts) { // merge v into u
                if(v>=u) continue; // parent should always have higher index
                auto nextSol = curSol;
                nextSol.merge(v,u);
                if(nextSol.width >= upper_bound) continue; // some pruning

                auto nextPar = curPar;
                rep(i,n) if(nextPar[i]==v) nextPar[i] = u;

                if(auto it = cache.find(nextPar); it==end(cache)) { // new reachable partition
                    cache[nextPar] = nextSol;
                    q.push(nextPar);
                } else if(nextSol.width < it->second.width) { // found better way to same partition
                    cache[nextPar] = nextSol;
                }
            }
        }
    }
    return {};
}

struct EdgeList {
    int n;
    vector<pair<int,int>> edges;
    vector<int> idx;
};

auto readInput() {
    string s = "c";
    while(s[0]!='p') getline(cin,s);
    int n,m;
    stringstream ss(s);
    ss >> s >> s >> n >> m;
    assert(s=="tww");
    vector<pair<int,int>> edges(m);
    for(auto& [u,v] : edges) {
        cin>>u>>v; --u; --v;
        assert(u!=v); // selfloop
    }
    vector<int> idx(n);
    iota(all(idx),0);
    return EdgeList{n,edges,idx};
}

auto components(const EdgeList& g) {
    // create adj list
    auto n = g.n;
    vector<vector<int>> adj(n);
    for(auto [u,v] : g.edges) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    // BFSs
    vector<int> comp(n,-1);
    int num_comps = 0;
    rep(i,n) {
        if(comp[i]!=-1) continue;
        int cur = num_comps++;
        comp[i] = cur;
        queue<int> q{{i}};
        while(size(q)) {
            int u = q.front(); q.pop();
            for(auto v : adj[u]) {
                if(comp[v]!=-1) continue;
                comp[v] = cur;
                q.push(v);
            }
        }
    }

    // create new edgelists
    vector<EdgeList> comps(num_comps);
    vector<int> local_idx(n); // index in own CC
    rep(i,n) local_idx[i] = comps[comp[i]].n++;
    rep(i,num_comps) comps[i].idx.resize(comps[i].n);
    rep(i,n) comps[comp[i]].idx[local_idx[i]] = i;
    for(auto [u,v] : g.edges)
        comps[comp[u]].edges.push_back({local_idx[u],local_idx[v]});
    return comps;
}

auto toMat(const EdgeList& g) {
    auto n = g.n;
    vector mat(n, vector(n, 0));
    for(auto [u,v] : g.edges) {
        mat[u][v] = mat[v][u] = 1;
    }
    return mat;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.precision(10);

    freopen("../../data/exact_050.gr", "r", stdin);

    auto edge_list = readInput();

    auto t1 = chrono::steady_clock::now();
    auto comps = components(edge_list);
    vector<Solution> sols;
    for(auto& comp : comps) {
        auto sol = partitionDP(toMat(comp));
        sols.push_back(*sol);
    }
    auto t2 = chrono::steady_clock::now();
    cout << "Branch and Bound Solution (";
    cout << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms)" << endl;
    int width = 0;
    for(auto& sol : sols) width = max(width,sol.width);
    cout << "twin width: " << width << endl;

    return 0;
}

