
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
    if(n<=20 || n>=3000) return 0;
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

int next_merge_bound(const Graph& mat) {
    if(size(mat)<=1 || size(mat)>=3000) return 0;
    int mn = 1e9;
    rep(i,ssize(mat)) rep(j,i) {
        if(mat[i][i]==-1 || mat[j][j]==-1) continue;
        mn = min(mn, reds_created(mat,j,i));
    }
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

int max_reds_after_merge(const Graph& mat, const vector<int>& rdeg, int u, int v) {
    assert(mat[u][u]!=-1);
    assert(mat[v][v]!=-1);
    auto n = ssize(mat);
    rep(i,n) assert(rdeg[i]==count(all(mat[i]),2));
    int own = 0;
    int mx_other = 0;
    rep(i,n) {
        if(i==u || i==v || mat[i][i]==-1) continue;
        auto [a,b] = minmax(mat[i][u], mat[i][v]);
        if(a == 2 || b==2 || pair(a,b) == pair(0,1)) own++;
        if(a == 2 && b==2) mx_other = max(mx_other, rdeg[i]-1);
        else if(a == 0 && b==1) mx_other = max(mx_other, rdeg[i]+1);
        else mx_other = max(mx_other, rdeg[i]);
    }
    return max(own, mx_other);
}

int next_bound2(const Graph& _mat) {
    if(size(_mat)>500) return 0;
    Graph mat = _mat;
    auto n = ssize(mat);
    vector<int> nodes;
    rep(i,n) if(mat[i][i]!=-1) nodes.push_back(i);
    vector<int> rdeg(n,0);
    rep(i,n) rdeg[i] = (int)count(all(mat[i]),2);

    int res = 0;
    while(size(nodes)>2) {
        // compute score
        sort(all(nodes));
        vector<vector<int>> score(n);
        for(int u : nodes)
            for(int v : nodes) {
                if(u>=v) continue;
                int badness = max_reds_after_merge(mat, rdeg, u, v);
                score[u].push_back(badness);
                score[v].push_back(badness);
            }
        for(auto& scores : score)
            sort(all(scores));

        // update bound
        sort(all(nodes), [&](int i, int j) { return score[i]>score[j]; });
        int next_bound = score[nodes.back()].front();
        res = max(res, next_bound);

        // del node with worst score
        int v = nodes.back(); nodes.pop_back();
        rdeg[v] = 0;
        rep(i,n) {
            if(mat[i][v]==2) rdeg[i]--;
            mat[i][v] = mat[v][i] = -1;
        }
        rep(i,n) assert(rdeg[i]==count(all(mat[i]),2));
//        cerr << next_bound << ' ';
    }
//    cerr << res << endl;
    return res;
}

int lower_bound(const Graph& mat, bool verbose) {
    auto g = aliveSubgraph(mat);

    auto t1 = chrono::steady_clock::now();
    int l0 = next_merge_bound(g);
    auto t2 = chrono::steady_clock::now();
    if(verbose) cerr << "next lb: " << l0 << " (" << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms)" << endl;

//    t1 = chrono::steady_clock::now();
//    int l1 = random_lower_bound(g, false);
//    t2 = chrono::steady_clock::now();
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

    t1 = chrono::steady_clock::now();
    int l6 = next_bound2(g);
    t2 = chrono::steady_clock::now();
    if(verbose) cerr << "next2lb: " << l6 << " (" << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << "ms)" << endl;

    // 91653349 // low = 0; no clusters; 10min
    // 23930461 // low = 2; no clusters; 2min
    // 11746941, 2732810 // low = 2; clusters; 126s
    //return max({l0,l1,l2,l3,l4,l5});

    // 34811736, 5321827
    // 35713037, 5427936

    // 23
    // 6033162; 88s; cluster
    // 6033162; 86s; no cluster
    return max({l0,l2,l3,l4,l5,l6});
}

