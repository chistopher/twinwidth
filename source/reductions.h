
#pragma once
#include "algo.h"

// 09 481687087
// 11 213779992 (80s)
// 15 823278744 (5min)
// 16 346687841
// 20 607000000
// 23

vector<pair<int,int>> cut_greed(const Solution& partial, int lb) {
    auto& mat=partial.mat;
    int n = int(size(mat));
    if(lb<=1) return {}; // does not work for tw1 or tw2 (e.g. path gets red at both ends)

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
            if(ssize(comp)+1<num_alive && ssize(comp)>2 && ssize(comp)<=3*lb) {
//                map<int,string> cols;
//                for(int vv : comp) cols[vv] = "blue";
//                cols[v] = "green";
//                draw(mat, cols);
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
                if(sol.width<=lb) {
                    vector<pair<int,int>> ret;
                    for(auto i=size(partial.merges.order); i<size(sol.merges.order); i++)
                        ret.emplace_back(sol.merges.order[i],sol.merges.parent[sol.merges.order[i]]);
                    return ret;
                }
            }
        }
    }
    return {};
}

vector<pair<int,int>> twins(const Solution& partial) {
    auto n = ssize(partial.mat);
    // find out nodes in graph merged until i-1
    vector<bool> alive(n,true);
    for(auto& x: partial.merges.order)
        alive[x] = false;

    vector used(n,false);
    vector<pair<int,int>> res;
    // merge twins
    rep(v,n) rep(u,n) { // v into u
        if(!alive[v] || !alive[u] || u==v) continue;
        if(used[v] || used[u]) continue;
        bool eq = true;
        rep(i,n) {
            if(!alive[i] || i==v || i==u) continue;
            if(minmax(partial.mat[v][i],partial.mat[u][i])==minmax(0,1)) eq = 0;
            if(partial.mat[v][i]==2 && partial.mat[u][i]!=2) eq = 0; // u has superset of red edges
            if(eq==0) break;
        }
        if(eq) {
            if(v>u) swap(v,u);
//            map<int,string> cols;
//            cols[v] = "red";
//            cols[u] = "red";
//            draw(partial.mat, cols);
            res.emplace_back(v,u);
            used[v] = used[u] = true;
            //return {{v,u}};
        }
    }
    return res;
}
