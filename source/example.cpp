
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

    struct MergeData {
        int width_before;
        vector<int> rdeg, adj_from, adj_to;
    };
    vector<MergeData> merge_data = {};

    void merge(int from, int to) {
        merge_data.push_back({width, rdeg, mat[from], mat[to]});
        rep(i,ssize(mat)) {
            if(i==from || i==to) continue;
            auto same = (mat[i][from] == mat[i][to]);
            auto new_val = same ? mat[i][to] : 2;
            //if(mat[i][to] == 2 && new_val != 2) rdeg[i]--, rdeg[to]--; // cannot happen
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

    void pop_merge() {
        auto n = ssize(mat);
        auto& prev = merge_data.back();
        int from = merges.order.back(); merges.order.pop_back();
        int to = merges.parent[from]; merges.parent[from] = -1;
        width = prev.width_before;
        rdeg = prev.rdeg;
        rep(i,n) {
            mat[i][to] = mat[to][i] = prev.adj_to[i];
            mat[i][from] = mat[from][i] = prev.adj_from[i];
        }
        merge_data.pop_back();
    }
};

struct Algo {
    bool verbose = false;
    const Graph g;
    map<vector<int>,int> cache;
    Solution best;
    int lower_bound = 0;
    int search_space = 0;

    explicit Algo(const Graph& _g) : g(_g) {};

    // runtime: fac(n) * fac(n-1) // 2**(n-1) * poly(n)
    void iterate_trees(Solution& partial, vector<int> partition) {
        if(best.width <= lower_bound) return; // we are done -> break out of recursion

        search_space++;
        if(verbose && search_space%100000 == 0) cout << "BnB search space: " << search_space << endl;

        auto n = ssize(partial.mat);
        if(ssize(partial.merges.order) == n - 1) { // complete tree
            if(partial.width < best.width) {
                best = partial;
                if(verbose) cout << "BnB found " << partial.width << "\t (after " << search_space << ')' << endl;
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
                iterate_trees(partial,partition);
                partial.pop_merge();
                return;
            }
        }

        vector<pair<int,int>> branches;

        rep(v,n) { // iterate possible candidates for order
            if(!alive[v]) continue;
            for(int p=v+1; p<n; ++p) { // iterate all possible parents
                if(!alive[p]) continue;
                branches.emplace_back(v,p);
            }
        }

        sort(all(branches),[&](auto& a, auto& b) { // merge red edges first :)
            return (partial.mat[a.first][a.second]==2) > (partial.mat[b.first][b.second]==2);
        });

        for(auto [v,p]: branches) {
            auto nextPar = partition;
            replace(all(nextPar),v,p);
            partial.merge(v,p);
            iterate_trees(partial,nextPar);
            partial.pop_merge();
        }
    }

    Solution solve(int lb = 0, int ub = 1e9) {
        // reset stuff
        cache.clear();
        best = Solution{};
        search_space = 0;

        auto n = ssize(g);
        best.width = ub;
        lower_bound = lb;
        Solution partial{g,{.order = {},.parent = vector(n,-1)}, 0, vector(n,0)};
        vector partition(n,0);
        iota(all(partition),0);
        iterate_trees(partial, partition);
        return best;
    }
};

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

Graph subgraph(const Graph& g, const vector<int>& nodes) {
    auto n = ssize(nodes);
    Graph mat(n, vector(n, 0));
    rep(i,n) rep(j,n) mat[i][j] = g[nodes[i]][nodes[j]];
    return mat;
}

int random_lower_bound(const Graph& mat, bool verbose) {
    if(size(mat)<=30) return 0;
    mt19937 gen{random_device{}()};
    vector perm(size(mat),0);
    iota(all(perm),0);
    int lower_bound = 0;
    rep(t,50000) {
        int N = 20;
        shuffle(all(perm),gen);
        auto nodes = perm;
        nodes.resize(N);
        auto g = subgraph(mat,nodes);
        auto sol = Algo(g).solve(lower_bound, lower_bound+2);
        if(sol.width > lower_bound) {
            if(verbose) cout << "found better lower bound: " << sol.width << endl;
            lower_bound = sol.width;
        }
    }
    return lower_bound;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.precision(10);

    freopen("../../data/exact_068.gr", "r", stdin);

    auto edge_list = readInput();
    cout << "n = " << edge_list.n << endl;
    cout << "deg = " << 2*ssize(edge_list.edges)/edge_list.n << endl;
    int l = 0;
    {
        auto t1 = chrono::steady_clock::now();
        l = random_lower_bound(toMat(edge_list), true);
        auto t2 = chrono::steady_clock::now();
        cout << "done computing lower bounds in " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms" << endl;
    }

    auto t1 = chrono::steady_clock::now();
    int search_space = 0;
    auto comps = components(edge_list);
    vector<Solution> sols;
    for(auto& comp : comps) {
        Algo alg(toMat(comp));
        alg.verbose = true;
        auto sol = alg.solve(l);
        sols.push_back(sol);
        search_space += alg.search_space;
    }
    auto t2 = chrono::steady_clock::now();
    cout << "Branch and Bound Solution" << endl;
    int width = 0;
    for(auto& sol : sols) width = max(width,sol.width);
    cout << "twin width:   " << width << endl;
    cout << "search space: " << search_space << endl;
    cout << "running time: " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms" << endl;

    return 0;
}

