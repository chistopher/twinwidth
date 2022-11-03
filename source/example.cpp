
#include <bits/stdc++.h>
#include <girgs/Generator.h>

#define rep(a, b)   for(int a = 0; a < (b); ++a)
#define all(a)      (a).begin(),(a).end()
//#define endl        '\n'

using namespace std;
using ll = long long;
using Graph = vector<vector<int>>; // adj matrix; 0: no edge 1: edge 2: red edge -1: deleted vertex

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
        rep(i,ssize(mat))
            mat[i][from] = mat[from][i] = -1;
        rep(i,ssize(mat))
            width = max(width, (int)count(all(mat[i]),2));
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

map<vector<int>,int> cache;
// runtime: fac(n) * fac(n-1) // 2**(n-1) * poly(n)
void iterate_trees(Solution partial, Solution& best, vector<int> partition) {
    auto n = ssize(partial.mat);
    if(size(partial.merges.order) == n - 1) { // complete tree
        if(partial.width < best.width) {
            best = partial;
            cout << "BnB found " << partial.width << endl;
        }
        return;
    }

    // early pruning
    if(partial.width>=best.width) return;
    if(auto it = cache.find(partition); it!=end(cache) && it->second<=partial.width)
        return;
    else cache[partition] = partial.width;

    // find out nodes in graph merged until i-1
    vector<bool> alive(n,true);
    for(auto& x: partial.merges.order)
        alive[x] = false;

    rep(v,n) { // iterate possible candidates for order
        if(!alive[v]) continue;
        for(int p=v+1; p<n; ++p) { // iterate all possible parents
            if(!alive[p]) continue;
            auto nxt = partial; nxt.merge(v,p);
            auto nextPar = partition;
            rep(i,n) if(nextPar[i]==v) nextPar[i] = p;
            iterate_trees(nxt,best,nextPar);
        }
    }
}

Solution branch_and_bound(const Graph& g) {
    auto n = ssize(g);
    cache.clear();
    Solution best, partial{g,{.order = {},.parent = vector(n,-1)},0};
    vector partition(n,0);
    iota(all(partition),0);
    iterate_trees(partial, best, partition);
    return best;
}

Graph generateER(int n, double p) {
    uniform_real_distribution dist;
    mt19937 gen;
    vector mat(n, vector(n, 0));
    rep(i,n) rep(j,n)
        if(i<j && dist(gen)<p)
            mat[i][j] = mat[j][i] = 1;
    return mat;
}


int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.precision(10);

    int n = 30;
    int d = 2;
    auto alpha = numeric_limits<double>::infinity();
    auto weights = girgs::generateWeights(n,2.8,213122412);
    auto positions = girgs::generatePositions(n,d,2242234);
    girgs::scaleWeights(weights,9,d,alpha);
    auto edge_list = girgs::generateEdges(weights,positions,alpha,2242-513);

    // print edge list
    for(auto& e: edge_list)
        cout << e.first << " " << e.second << endl;
    // make mat matrix
    vector mat(n, vector(n, 0));
    for(auto& [u,v]: edge_list)
        mat[u][v] = mat[v][u] = 1;

    mat = generateER(n,0.3);

    auto t1 = chrono::steady_clock::now();
    auto sol = branch_and_bound(mat);
    auto t2 = chrono::steady_clock::now();
    cout << "Branch and Bound Solution (";
    cout << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms)" << endl;
    cout << "twin width: " << sol.width << endl;
    cout << "merges: " << endl;
    for(auto& x: sol.merges.order)
        cout << '\t' << x << " into " << sol.merges.parent[x] << endl;

    return 0;
}

