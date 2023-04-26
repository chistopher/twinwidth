
#include "lower_bounds.h"

#include "algo.h"
#include "graph.h"

int random_lower_bound(const Graph& mat, bool verbose) {
    if(size(mat)<=30) return 0;
    //mt19937 gen{random_device{}()};
    mt19937 gen{378};
    vector perm(size(mat),0);
    iota(all(perm),0);
    int lower_bound = 0;
    rep(t,50000) {
        int N = 20;
        shuffle(all(perm),gen);
        auto nodes = perm;
        nodes.resize(N);
        auto g = subgraph(mat,nodes);
        auto sol = Algo().solve(g,lower_bound, lower_bound+2);
        if(sol.width > lower_bound) {
            if(verbose) cerr << "found better lower bound: " << sol.width << endl;
            lower_bound = sol.width;
        }
    }
    return lower_bound;
}

// thomas page rank lower bound
int lower_bound2(const Graph& mat) {
    int n = (int)size(mat);
    if(n>200) return 0;
    vector<double> rank(n,1.0), next(n);
    rep(t,100) {
        fill(all(next),0.0);
        for(int i=0; i<n; ++i)
            for(int j=i+1; j<n; ++j)
                rep(k,n)
                    if(k!=i && k!=j && mat[i][k]!=mat[j][k])
                        next[k] += (rank[i] * rank[j]);
        auto sum = accumulate(all(next),0.0);
        for(auto& x : next) x /= sum;
        swap(rank,next);
    }
    vector<int> nodes;
    { // find top 25 connected nodes
        set<pair<double,int>> q;
        int s = int(max_element(all(rank)) - begin(rank));
        q.emplace(rank[s],s);
        vector done(n,0);
        done[s] = 1;
        while(size(q)) {
            auto [r,u] = *q.rbegin(); q.erase(prev(end(q)));
            nodes.push_back(u);
            for(int v=0; v<n; ++v) {
                if(!done[v] && mat[u][v]>0)
                    q.emplace(rank[v],v), done[u]=1;
            }
        }
        nodes.resize(min(n,25));
    }
    auto l = Algo().solve(subgraph(mat,nodes)).width;
    return l;
}

int bunk_bound(const Graph& mat) {
    auto n = ssize(mat);
    if(n<=20) return 0;
    int low = 0;
    mt19937 gen{};
    rep(s,n) {
        vector<int> q{s};
        vector<int> dist(n,-1);
        dist[s] = 0;
        auto done = 0u;
        while(done<size(q)) {
            int v = q[done++];
            rep(u,n) {
                if(dist[u]!=-1 || mat[v][u]==0) continue;
                dist[u] = dist[v]+1;
                q.push_back(u);
            }
        }

        // permute nodes inside a BFS layer
        unsigned l=0, r=0;
        while(r<done) {
            while(r<done && dist[q[l]]==dist[q[r]]) r++;
            shuffle(q.begin()+l, q.begin()+r, gen);
            l = r;
        }

        // brute force subgraph around node s
        auto nodes = q;
        nodes.resize(20);
        low = max(low, Algo().solve(subgraph(mat,nodes)).width);
    }
    return low;
}

int lower_bound(const Graph& mat, bool verbose) {
    auto g = aliveSubgraph(mat);
    auto t1 = chrono::steady_clock::now();
    int l1 = random_lower_bound(g, false);
    auto t2 = chrono::steady_clock::now();
    if(verbose) cerr << "rand lb: " << l1 << " (" << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms)" << endl;
    t1 = chrono::steady_clock::now();
    auto l2 = lower_bound2(g);
    t2 = chrono::steady_clock::now();
    if(verbose) cerr << "page lb: " << l2 << " (" << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms)" << endl;
    t1 = chrono::steady_clock::now();
    auto l3 = bunk_bound(g);
    t2 = chrono::steady_clock::now();
    if(verbose) cerr << "bunk lb: " << l3 << " (" << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms)" << endl;
    return max({l1,l2,l3});
}

