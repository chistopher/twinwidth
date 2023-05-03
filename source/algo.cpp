
#include "algo.h"
#include "reductions.h"
#include "lower_bounds.h"

void Solution::merge(int from, int to) {
    assert(from<to);
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

void Solution::pop_merge() {
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

void Algo::iterate_trees(Solution& partial, vector<int> partition) {
    if(best.width <= m_lower_bound) return; // we are done -> break out of recursion

    search_space++;
    static auto min_seen = size(partial.merges.order);
    min_seen = min(min_seen, size(partial.merges.order));
    if(verbose && search_space%10000 == 0) cerr << "BnB search space: " << search_space << ' ' << min_seen << endl, min_seen = 10000;

    auto n = ssize(partial.mat);
    if(ssize(partial.merges.order) == n - 1) { // complete tree
        if(partial.width < best.width) {
            best = partial;
            if(verbose) cerr << "BnB found " << partial.width << "\t (after " << search_space << ')' << endl;
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
            if(v>u) swap(v,u);
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

    vector<int> badness(size(branches),0);
    vector<int> rdeg(size(branches),0);
    vector<int> nodes;
    rep(i,n) if(alive[i]) nodes.push_back(i);
    int INF = 1000000;
    rep(t,ssize(branches)) {
        auto [u,v] = branches[t];
        auto& mat = partial.mat;
        if(mat[u][v]==2) badness[t]--;
        for(int i : nodes) {
            if(i==u || i==v) continue;
            if(mat[i][u]==2 || mat[i][v]==2 || minmax(mat[i][u],mat[i][v])==minmax(0,1)) rdeg[t]++;
            if(minmax(mat[u][i],mat[v][i])==minmax(0,1)) badness[t]++;
            else if(mat[u][i]==2 && mat[v][i]==2) badness[t]--;
            if(rdeg[t]>=best.width) { badness[t] = INF; break; }
        }
    }
    vector<int> perm(ssize(branches));
    iota(all(perm),0);
    sort(all(perm),[&](auto& a, auto& b) { // merge red edges first :)
        return badness[a] < badness[b];
        //return (partial.mat[a.first][a.second]==2) > (partial.mat[b.first][b.second]==2);
    });
    while(!empty(perm) && badness[perm.back()]==INF) perm.pop_back();

    rep(i,ssize(perm)) {
        if(rdeg[perm[i]]>=best.width) continue;
        auto [v,p] = branches[perm[i]];
        auto nextPar = partition;
        replace(all(nextPar),v,p);
        partial.merge(v,p);
        iterate_trees(partial,nextPar);
        partial.pop_merge();
        if(partial.width>=best.width) return;
    }
}

Solution Algo::solve(const Graph& g, int lb, int ub) {
    rep(i,ssize(g)) assert(g[i][i]!=-1);
    // reset stuff
    cache.clear();
    best = Solution{};
    search_space = 0;

    auto n = ssize(g);
    best.width = ub;
    m_lower_bound = lb;
    Solution partial{g,{.order = {},.parent = vector(n,-1)}, lb, vector(n,0)};
    rep(i,n) partial.rdeg[i] = (int)count(all(g[i]),2); // needed if called with red edges

    if(verbose) {
        vector<pair<int,int>> m;
        do {
            string rule;
            m = twins(partial);
            if(!empty(m)) rule = "twin";
            if(empty(m)) m=cut_greed(partial,lb);
            if(!empty(m) && empty(rule)) rule = "greed" + to_string(size(m));
            for(auto [v,p] : m)
                partial.merge(v,p);
            if(size(m)) cerr << "reduced to " << ssize(partial.mat) - ssize(partial.merges.order) << " with " << rule << endl;
        } while(!empty(m));
    }

    if(partial.merges.order.size()) { // initial reductions worked
        auto l2 = lower_bound(partial.mat, true);
        m_lower_bound = lb = max(l2, lb);
        // remove deleted nodes
        vector<int> nodes;
        rep(i,n) if(partial.mat[i][i]!=-1) nodes.push_back(i);
        auto g2 = subgraph(partial.mat, nodes);
        auto rdeg = vector(ssize(g2),0);
        rep(i,ssize(g2)) rep(j,i) if(g2[i][j]==2) rdeg[i]++, rdeg[j]++;
        auto partial2 = Solution{g2,{.order = {},.parent = vector(ssize(g2),-1)}, lb, rdeg};
        vector partition2(size(g2),0);
        iota(all(partition2),0);
        iterate_trees(partial2, partition2);
        if(best.width>=1'000'000'000) return Solution{};
        // apply initial merges to solution
        rep(i,ssize(best.merges.order)) {
            int v = best.merges.order[i];
            int p = best.merges.parent[v];
            partial.merge(nodes[v],nodes[p]);
        }
        return partial;
    }

    // no initial reductions
    vector partition(n,0);
    iota(all(partition),0);
    for(int v : partial.merges.order)
        replace(all(partition),v,partial.merges.parent[v]);
    iterate_trees(partial, partition);
    return best;
}
