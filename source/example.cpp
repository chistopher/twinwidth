
#include <bits/stdc++.h>

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
    vector<int> rdeg;
    void merge(int from, int to) {
        rep(i,ssize(mat)) {
            if(i==from || i==to) continue;
            auto same = (mat[i][from] == mat[i][to]);
            auto new_val = same ? mat[i][to] : 2;
            if(mat[i][to] == 2 && new_val != 2) rdeg[i]--, rdeg[to]--;
            if(mat[i][to] != 2 && new_val == 2) rdeg[i]++, rdeg[to]++;
            mat[i][to] = mat[to][i] = new_val;

        }
        rep(i,ssize(mat)) {
            if(mat[i][from] == 2) rdeg[i]--, rdeg[from]--;
            mat[i][from] = mat[from][i] = -1;
        }
        //rep(i,ssize(mat)) width = max(width, (int)count(all(mat[i]),2));
        width = max(width, *max_element(all(rdeg)));
        rep(i,ssize(mat)) assert(rdeg[i]==count(all(mat[i]),2));
        merges.order.push_back(from);
        merges.parent[from] = to;
    }
};

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

    // merge twins
    rep(v,n) rep(u,n) { // v into u
        if(!alive[v] || !alive[u] || u==v) continue;
        bool eq = 1;
        //rep(i,n) if(alive[i] && i!=v && i!=u && partial.mat[v][i]!=partial.mat[u][i]) eq = 0;
        rep(i,n) {
            if(!alive[i] || i==v || i==u) continue;
            if(minmax(partial.mat[v][i],partial.mat[u][i])==minmax(0,1)) eq = 0;
            if(partial.mat[v][i]==2 && partial.mat[u][i]!=2) eq = 0;
            if(eq==0) break;
        }
        if(eq) {
            partial.merge(v,u);
            replace(all(partition),v,u);
            return iterate_trees(move(partial),best,partition);
        }
    }

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
    Solution best, partial{g,{.order = {},.parent = vector(n,-1)},0, vector(n,0)};
    vector partition(n,0);
    iota(all(partition),0);
    iterate_trees(partial, best, partition);
    return best;
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

    //freopen("../../data/exact_050.gr", "r", stdin);

    auto edge_list = readInput();

    auto t1 = chrono::steady_clock::now();
    auto comps = components(edge_list);
    vector<Solution> sols;
    for(auto& comp : comps) {
        auto sol = branch_and_bound(toMat(comp));
        sols.push_back(sol);
    }
    auto t2 = chrono::steady_clock::now();
    cout << "Branch and Bound Solution (";
    cout << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms)" << endl;
    int width = 0;
    for(auto& sol : sols) width = max(width,sol.width);
    cout << "twin width: " << width << endl;

    return 0;
}

