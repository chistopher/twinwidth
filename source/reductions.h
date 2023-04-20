
#pragma once
#include "algo.h"

// 09 481687087
// 11 213779992
// 16 346687841

optional<Solution> cut_greed(const Solution& partial, int lb) {
    auto& mat=partial.mat;
    int n = int(size(mat));

    // find out nodes in graph merged until i-1
    vector<bool> alive(n,true);
    for(auto& x: partial.merges.order)
        alive[x] = false;

    int num_alive = accumulate(all(alive),0);
    rep(v,n) {
        if(!alive[v]) continue;
        vector dist(n,-1);
        dist[v] = 0;
        rep(nei,n) {
            if(!alive[nei] || dist[nei]!=-1 || mat[v][nei]==0) continue;
            // start bfs from mat
            vector<int> comp{nei};
            dist[nei] = 1;
            auto done = 0u;
            while(done<size(comp)) {
                int u = comp[done++];
                rep(w,n) {
                    if(!alive[w] || dist[w]!=-1 || mat[u][w]==0) continue;
                    dist[w] = dist[u]+1;
                    comp.push_back(w);
                }
            }
            if(ssize(comp)+1<num_alive && ssize(comp)>2 && ssize(comp)<=lb) {
                Solution sol = partial;
                vector<int> dist1, other;
                for(auto u : comp) if(mat[v][u]!=0) dist1.push_back(u); else other.push_back(u);
                if(!empty(other)) {
                    int mx_other = *max_element(all(other));
                    for(auto u : other)
                        if(u!=mx_other)
                            sol.merge(u,mx_other);
                }
                int mx_dist1 = *max_element(all(dist1));
                for(auto u : dist1)
                    if(u!=mx_dist1)
                        sol.merge(u,mx_dist1);
                assert(sol.merges.order>partial.merges.order);
                if(sol.width<=lb) return sol;
                //cout << "found comp of size " << ssize(comp) << " and lb " << lb << endl;
            }
        }
    }
    return {};
}
