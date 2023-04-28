
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

// other page rank lb
int lower_bound3(const Graph& mat) {
    int n = (int)size(mat);
    if(n>200) return 0;
    vector<double> rank(n,1.0), next(n);
    rep(t,100) {
        fill(all(next),0.0);
        rep(i,n) rep(j,i) {
            double diff = 0;
            rep(k,n) if(k!=i && k!=j && mat[i][k]!=mat[j][k]) diff += rank[k];
            next[i] += diff;
            next[j] += diff;
        }
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

int reds_created(const Graph& mat, int u, int v) {
    assert(mat[u][u]!=-1);
    assert(mat[v][v]!=-1);
    auto n = ssize(mat);
    int res = 0;
    rep(i,n) {
        if(i==u || i==v || mat[i][i]==-1) continue;
        if(minmax(mat[u][i],mat[v][i])==minmax(0,1)) res++;
        else if(mat[u][i]==2 && mat[v][i]==2) res--;
    }
    return res - (mat[u][v]==2);
}

int foo(const Graph& mat) {
    if(size(mat)<=1) return 0;
    int mn = 1e9;
    rep(i,ssize(mat)) rep(j,i) mn = min(mn, reds_created(mat,j,i));
    return mn;
}

int lower_bound4(const Graph& mat) {
    auto n = ssize(mat);
    int res = 0;
    int max_v = 0;
    rep(i,n) if(count(all(mat[i]),1)>count(all(mat[max_v]),1)) max_v = i;
    rep(t,n) {
        if(t!=max_v) continue;
        // seed
        vector nodes{t};
        rep(i,n) if(i!=t && mat[t][i]>=1 && size(nodes)<5) nodes.push_back(i);

        vector inside(n,false);
        for(auto x : nodes) inside[x] = true;

        // grow
        while(size(nodes)<25) {
            int best = -1;
            int best_score = 0;
            rep(i,n) {
                if(inside[i]) continue;
                int score = 1e9;
                for(auto x : nodes) {
                    int merge_cost = 0; // simulated merge of x and i
                    for(auto y : nodes) merge_cost += (mat[i][y]!=mat[x][y]);
                    score = min(score,merge_cost);
                    if(score<=best_score) break;
                }
                if(score>best_score) best_score=score, best=i;
            }
            if(best==-1) break;
            nodes.push_back(best);
            inside[best] = true;
        }

        // test
        res = max(res,Algo().solve(subgraph(mat,nodes)).width);
    }
    return res;
}

int lower_bound(const Graph& mat, bool verbose) {
    auto g = aliveSubgraph(mat);
    int l0 = foo(g);
    if(verbose) cerr << "stunk b: " << l0 << endl;
    auto t1 = chrono::steady_clock::now();
//    int l1 = random_lower_bound(g, false);
    auto t2 = chrono::steady_clock::now();
//    if(verbose) cerr << "rand lb: " << l1 << " (" << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms)" << endl;
    t1 = chrono::steady_clock::now();
    auto l2 = lower_bound2(g);
    t2 = chrono::steady_clock::now();
    if(verbose) cerr << "page lb: " << l2 << " (" << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms)" << endl;
    t1 = chrono::steady_clock::now();
    auto l3 = bunk_bound(g);
    t2 = chrono::steady_clock::now();
    if(verbose) cerr << "bunk lb: " << l3 << " (" << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms)" << endl;

    t1 = chrono::steady_clock::now();
    auto l4 = lower_bound3(g);
    t2 = chrono::steady_clock::now();
    if(verbose) cerr << "lb3:     " << l4 << " (" << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms)" << endl;

    t1 = chrono::steady_clock::now();
    auto l5 = lower_bound4(g);
    t2 = chrono::steady_clock::now();
    if(verbose) cerr << "lb4:     " << l5 << " (" << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms)" << endl;

    //return max({l0,l1,l2,l3,l4,l5});
    return max({l0,l2,l3,l4,l5});
}

