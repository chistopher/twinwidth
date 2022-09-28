
#include <bits/stdc++.h>
#include <girgs/Generator.h>

#define rep(a, b)   for(int a = 0; a < (b); ++a)
#define all(a)      (a).begin(),(a).end()
#define endl        '\n'

using namespace std;
using Graph = vector<vector<int>>;
using ll = long long;

struct TreeRepresentation {
    vector<int> order;
    vector<int> parent;
};


// runtime: fac(n) * fac(n-1) // 2**(n-1)
template<typename T>
void iterate_trees(int n, TreeRepresentation& partial, T&& callback) {
    if(size(partial.order)==n-1) { // complete tree
        callback(partial);
        return;
    }

    // early pruning
    // reduction rules

    // find out nodes in graph merged until i-1
    vector<bool> alive(n,true);
    for(auto& x: partial.order)
        alive[x] = false;

    rep(v,n) { // iterate possible candidates for order
        if(!alive[v]) continue;
        partial.order.push_back(v);
        for(int p=v+1; p<n; ++p) { // iterate all possible parents
            if(!alive[p]) continue;
            partial.parent[v] = p;
            iterate_trees(n,partial,callback);
        }
        partial.order.pop_back();
    }

}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.precision(10);

    int n = 9;
    int d = 2;
    auto alpha = numeric_limits<double>::infinity();
    auto weights = girgs::generateWeights(n,2.8,213122412);
    auto positions = girgs::generatePositions(n,d,2242234);
    girgs::scaleWeights(weights,4,d,alpha);
    auto edge_list = girgs::generateEdges(weights,positions,alpha,2242-513);

    // print edge list
    for(auto& e: edge_list)
        cout << e.first << " " << e.second << endl;
    // make mat matrix
    vector mat(n, vector(n, 0));
    for(auto& [u,v]: edge_list)
        mat[u][v] = mat[v][u] = 1;

    int counter = 0;
    int twin_width = n+1, num_found = 0;
    TreeRepresentation solution;

    TreeRepresentation partial;
    partial.parent.resize(n);
    iterate_trees(n,partial,[&](auto& tree) {
        counter++;
        auto graph = mat;
        vector<bool> alive(n,true);

        auto max_deg = 0;

        // do merges
        rep(i,n-1) {
            auto v = tree.order[i];
            auto p = tree.parent[v];
            assert(alive[v] && alive[p]);

            rep(u,n) {
                if(!alive[u]) continue;
                auto same = (graph[v][u] == graph[p][u]);
                auto new_val = same ? graph[v][u] : 2;
                graph[p][u] = graph[u][p] = new_val;
            }

            rep(u,n)
                if(alive[u])
                    max_deg = max(max_deg, (int)count(all(graph[u]),2));

            alive[v] = false;
            rep(u,n) graph[v][u] = graph[u][v] = -1;
        }


        if(twin_width > max_deg) {
            cout << counter << " " << max_deg << ' ' << num_found << endl;
            cout.flush();
            twin_width = max_deg;
            solution = tree;
            num_found = 1;
        } else if(twin_width == max_deg) {
            num_found++;
        }
    });
    cout << counter << ' ' << num_found << endl;
    cout << "twin width: " << twin_width << endl;
    cout << "order: " << endl;
    for(auto& x: solution.order)
        cout << x << " ";
    cout << n-1 << endl;
    cout << "parents: " << endl;
    for(auto& x: solution.parent)
        cout << x << " ";
    cout << endl;



    return 0;
}

