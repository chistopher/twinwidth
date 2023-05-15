
#include "algo.h"
#include "reductions.h"
#include "lower_bounds.h"

vector<int> Solution::TreeRepresentation::get_partition() const {
    auto n = ssize(parent);
    auto partition = parent;
    for(auto i=n-1; i>=0; --i)
        partition[i] = int(parent[i]==-1 ? i : partition[parent[i]]);
    return partition;
}

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
        //draw(partial.mat);
        if(verbose) cerr << "reductions done" << endl;
    }

    if(empty(partial.merges.order)) { // no initial reductions
        iterate_trees(partial);
        return best;
    }

    // initial reductions worked
    auto l2 = lower_bound(partial.mat, true);
    m_lower_bound = lb = max(l2, lb);
    // remove deleted nodes
    vector<int> nodes;
    rep(i,n) if(partial.mat[i][i]!=-1) nodes.push_back(i);
    auto g2 = subgraph(partial.mat, nodes);
    auto rdeg = vector(ssize(g2),0);
    rep(i,ssize(g2)) rep(j,i) if(g2[i][j]==2) rdeg[i]++, rdeg[j]++;
    auto partial2 = Solution{g2,{.order = {},.parent = vector(ssize(g2),-1)}, lb, rdeg};
    iterate_trees(partial2); // solve subgraph after reductions
    // now solution to reduced instance is in `best` and we need to lift it to a solution for the original instance
    if(best.width>=1'000'000'000) return Solution{};
    rep(i,ssize(best.merges.order)) {
        int v = best.merges.order[i];
        int p = best.merges.parent[v];
        partial.merge(nodes[v],nodes[p]);
    }
    return partial;
}

void Algo::iterate_trees(Solution& partial) {
    if(best.width <= m_lower_bound) return; // we are done -> break out of recursion

    search_space++;
    static auto min_seen = size(partial.merges.order);
    min_seen = min(min_seen, size(partial.merges.order));
    if(verbose && search_space%10000 == 0) cerr << "BnB search space: " << search_space << ' ' << min_seen << endl, min_seen = 10000;

    auto n = ssize(partial.mat);
    if(ssize(partial.merges.order) == n - 1) { // complete tree
        if(partial.width < best.width) {
            best = partial;
            cluster.clear();
            if(verbose) cerr << "BnB found " << partial.width << "\t (after " << search_space << ')' << endl;
//            if(verbose) {
//                auto sol = best;
//                while(sol.merges.order.size()) {
//                    int v = sol.merges.order.back();
//                    int p = sol.merges.parent[v];
//                    sol.pop_merge();
//                    draw_snapshot(sol.mat, sol.merges.order.size(), v, p);
//                }
//            }
        }
        return;
    }

    // early pruning
    const auto partition = partial.merges.get_partition();
    if(partial.width>=best.width) return;
    if(auto it = cache.find(partition); it!=end(cache) && it->second<=partial.width)
        return;
    else cache[partition] = partial.width;

//    if(verbose && !empty(partial.merges.order)) {
//        add_cluster(partial, partial.merges.parent[partial.merges.order.back()]);
//    }

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
            iterate_trees(partial);
            partial.pop_merge();
            return;
        }
    }

//    auto reduction = check_clusters(partial);
//    if(!empty(reduction)) {
//        //if(verbose) cerr << "reduction by clusters" << endl;
//        bool skip = false;
//        bool first = true;
//        for(auto [v,p] : reduction) {
//            if(!first) {
//                const auto partition = partial.merges.get_partition();
//                if(auto it = cache.find(partition); it!=end(cache) && it->second<=partial.width)
//                    skip = true;
//                else cache[partition] = partial.width;
//            }
//            first = false;
//            partial.merge(v,p);
//        }
//        if(!skip) iterate_trees(partial);
//        for(auto _ : reduction) partial.pop_merge();
//        return;
//    }

    // find and sort all merges
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

    // branch on all possible merges
    rep(i,ssize(perm)) {
        if(rdeg[perm[i]]>=best.width) continue;
        auto [v,p] = branches[perm[i]];
        partial.merge(v,p);
        iterate_trees(partial);
        partial.pop_merge();
        if(partial.width>=best.width) return;
    }
}

void Algo::add_cluster(const Solution &partial, int p) {
    auto n = ssize(partial.mat);
    auto partition = partial.merges.get_partition();

    vector sz(n,0);
    rep(i,n) sz[partition[i]]++;
    int red_con = 0;
    rep(i,n) if(partial.mat[p][i]==2) red_con += sz[i];
    if(red_con<=2 && sz[p]>4) {
        vector<int> nodes;
        rep(i,n) if(partition[i]==p) nodes.push_back(i);
        vector<pair<int,int>> merges;
        for(int x : partial.merges.order) {
            auto y = partial.merges.parent[x];
            if(partition[x]==p && partition[y]==p)
                merges.push_back({x,y});
        }
        assert(merges.size()==size(nodes)-1);
        cluster.insert({nodes, merges});

//        auto tmp = partial;
//        while(tmp.merges.order.size()) tmp.pop_merge();
//        map<int,string> colors;
//        rep(i,n) {
//            if(partition[i]==p) continue;
//            int rep = i;
//            while(partial.merges.parent[rep]!=-1) rep = partial.merges.parent[rep];
//            if(partial.mat[rep][p]==2) colors[i] = "red";
//        }
//        for(auto x : nodes) colors[x] = "blue";
//        draw(tmp.mat, colors);
    }
}

vector<pair<int, int>> Algo::check_clusters(Solution &partial) {
    int allowed_width = max(partial.width,m_lower_bound);
    auto n = ssize(partial.mat);
    auto super = [&](int x) {
        while(partial.merges.parent[x]!=-1) x = partial.merges.parent[x];
        return x;
    };

    for(auto [nodes, merges] : cluster) {

        bool worked = true;
        vector<pair<int,int>> applied_merges;
        for(auto [x,y] : merges) {
            int rep1 = super(x), rep2 = super(y);
            if(rep1==rep2) continue;
            if(rep1>rep2) swap(rep1,rep2);
            partial.merge(rep1,rep2);
            applied_merges.emplace_back(rep1,rep2);
            worked &= (partial.width<=allowed_width);
        }
        if(!worked || empty(applied_merges)) {
            for(auto _ : applied_merges) partial.pop_merge();
            continue;
        }

        auto cluster_rep = super(nodes.back());

        // check tree like stuff
        if(partial.rdeg[cluster_rep]==1 && allowed_width>=2) {
            int nei = -1;
            rep(i,n) if(partial.mat[cluster_rep][i]==2) nei = i;
            assert(nei!=-1);
            auto nei_deg_b = count(all(partial.mat[nei]),1);;
            auto nei_deg_r = count(all(partial.mat[nei]),2);
            if((nei_deg_r==1 && nei_deg_b==1) || (nei_deg_r==2 && nei_deg_b==0)) {
                for(auto _ : applied_merges) partial.pop_merge();
                return applied_merges;
            }
        }

        // check if there is a cluster node before merge that has red edges to entire red neighborhood of cluster after merge
        vector has_red_after_merge(n,false);
        rep(i,n) if(partial.mat[cluster_rep][i]==2) has_red_after_merge[i] = true;
        for(auto _ : applied_merges) partial.pop_merge();

        int has_all_node = -1;
        for(int x : nodes) {
            x = super(x);
            bool is_good = true;
            rep(i,n) if(has_red_after_merge[i]) is_good &= (partial.mat[x][i]==2);
            if(is_good) has_all_node = x;
        }
        if(has_all_node!=-1)
            return applied_merges; // safe
    }

    return {};
}
