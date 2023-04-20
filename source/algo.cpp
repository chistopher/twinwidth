
#include "algo.h"
#include "reductions.h"

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
    if(best.width <= lower_bound) return; // we are done -> break out of recursion

    search_space++;
    if(verbose && search_space%1000000 == 0) cerr << "BnB search space: " << search_space << endl;

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

    if(verbose && false) {
        auto s = cut_greed(partial,max(lower_bound,partial.width));
        if(s) {
            auto sol = *s;
            vector part = partition;
            for(auto& x: sol.merges.order)
                alive[x] = false;
            rep(i,n) while(!alive[part[i]]) part[i] = sol.merges.parent[part[i]];
            return iterate_trees(sol,part);
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

Solution Algo::solve(int lb, int ub) {
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
