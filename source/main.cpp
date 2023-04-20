
#include "graph.h"
#include "algo.h"
#include "lower_bounds.h"

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.precision(10);

    //freopen("../../data/exact_016.gr", "r", stdin);

    auto edge_list = readInput();
    cerr << "n = " << edge_list.n << endl;
    cerr << "deg = " << 2*ssize(edge_list.edges)/edge_list.n << endl;

    auto t1 = chrono::steady_clock::now();
    int search_space = 0;
    auto comps = components(edge_list);
    vector<Solution> sols;
    int l_max = 0;
    for(auto& comp : comps) {
        while(comp.n>=200); // TLE on large instances
        if(comps.size()>1) cerr << "component " << &comp - &comps[0] << endl;
        int l = lower_bound(toMat(comp), true);
        cerr << "lower bound:  " << l << endl;
        Algo alg(toMat(comp));
        alg.verbose = true;
        auto sol = alg.solve(l);
        sols.push_back(sol);
        search_space += alg.search_space;
        l_max = max(l_max, l);
    }

    auto t2 = chrono::steady_clock::now();
    cerr << "Branch and Bound Solution" << endl;
    int width = 0;
    for(auto& sol : sols) width = max(width,sol.width);
    cerr << "twin width:   " << width << endl;
    cerr << "search space: " << search_space << endl;
    cerr << "running time: " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " ms" << endl;

    // print output sequence
    rep(i,ssize(comps))
        for(int v : sols[i].merges.order)
            cout << comps[i].idx[v]+1 << ' ' << comps[i].idx[sols[i].merges.parent[v]]+1 << endl;

        // only validate in debug mode
#ifndef NDEBUG
    Solution sol{toMat(edge_list),{.order = {},.parent = vector(edge_list.n,-1)}, 0, vector(edge_list.n,0)};
    rep(i,size(comps))
        for(int v : sols[i].merges.order)
            sol.merge(comps[i].idx[v], comps[i].idx[sols[i].merges.parent[v]]);
    assert(sol.width==width);
    assert(ssize(sol.merges.order)==edge_list.n-1);
#endif

    if(argc>1 && string(argv[1])=="-v") {
        ofstream file("run_info.json");
        file << "{\n";
        file << "\t\"width\": " << width << ",\n";
        file << "\t\"search_space\": " << search_space << ",\n";
        file << "\t\"running_time\": " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << ",\n";
        file << "\t\"lower_bound\": " << l_max << "\n";
        file << "}";
    }

    return 0;
}

