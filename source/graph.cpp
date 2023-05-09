
#include "graph.h"

EdgeList readInput() {
    string s = "c";
    while(s[0]!='p') getline(cin,s);
    int n,m;
    stringstream ss(s);
    ss >> s >> s >> n >> m;
    assert(s=="tww");
    vector<pair<int,int>> edges(m);
    for(auto& [u,v] : edges) {
        cin>>u>>v; --u; --v;
        assert(u!=v); // selfloop
    }
    vector<int> idx(n);
    iota(all(idx),0);
    return {n,edges,idx};
}

vector<EdgeList> components(const EdgeList& g) {
    // create adj list
    auto n = g.n;
    vector<vector<int>> adj(n);
    for(auto [u,v] : g.edges) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    // BFSs
    vector<int> comp(n,-1);
    int num_comps = 0;
    rep(i,n) {
        if(comp[i]!=-1) continue;
        int cur = num_comps++;
        comp[i] = cur;
        queue<int> q{{i}};
        while(size(q)) {
            int u = q.front(); q.pop();
            for(auto v : adj[u]) {
                if(comp[v]!=-1) continue;
                comp[v] = cur;
                q.push(v);
            }
        }
    }

    // create new edgelists
    vector<EdgeList> comps(num_comps);
    vector<int> local_idx(n); // index in own CC
    rep(i,n) local_idx[i] = comps[comp[i]].n++;
    rep(i,num_comps) comps[i].idx.resize(comps[i].n);
    rep(i,n) comps[comp[i]].idx[local_idx[i]] = i;
    for(auto [u,v] : g.edges)
        comps[comp[u]].edges.push_back({local_idx[u],local_idx[v]});
    return comps;
}

Graph toMat(const EdgeList& g) {
    auto n = g.n;
    vector mat(n, vector(n, 0));
    for(auto [u,v] : g.edges) {
        mat[u][v] = mat[v][u] = 1;
    }
    return mat;
}

Graph subgraph(const Graph& g, const vector<int>& nodes) {
    auto n = ssize(nodes);
    Graph mat(n, vector(n, 0));
    rep(i,n) rep(j,n) mat[i][j] = g[nodes[i]][nodes[j]];
    return mat;
}


Graph aliveSubgraph(const Graph& g) {
    vector<int> nodes;
    rep(i,ssize(g)) if(g[i][i]!=-1) nodes.push_back(i);
    return subgraph(g,nodes);
}


void draw(const Graph& g, map<int,string> color, map<int,string> label) {
    // use graphviz to draw graph
    ofstream f("graph.dot");
    f << "graph G {" << endl;
    rep(i,ssize(g)) {
        if(g[i][i]==-1) continue; // deleted vertex
        f << i;
        if(color.count(i)) f << " [style=filled,fillcolor=" << color[i] << "]";
        if(label.count(i)) f << " [label=" << label[i] << "]";
        f << ";\n" << endl;
    }
    rep(i,ssize(g)) rep(j,i) {
        if(g[i][j]==0 || g[i][j]==-1) continue;
        f << i << " -- " << j;
        if(g[i][j]==2) f << " [color=red]";
        f << ";\n" << endl;
    }
    f << "}" << endl;
    system("dot -Tpng graph.dot -o graph.png");
}

void draw_snapshot(const Graph& g, int time, int v, int p) {
    int n = size(g);
    std::ofstream file(to_string(time) + ".gexf");
    file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    file << "<gexf version=\"1.3\" xmlns=\"http://gexf.net/1.3\" xmlns:viz=\"http://gexf.net/1.3/viz\">\n";
    file << "<graph mode=\"slice\" defaultedgetype=\"undirected\" timerepresentation=\"timestamp\" timestamp=\"" << time << "\">\n";
    file << "<nodes>\n";
    for(int i = 0; i < n; ++i)
        if(g[i][i]!=-1) {
            file << "\t<node id=\"" << i << "\">";
            if(i==v) file << "<viz:color r=\"9\" g=\"200\" b=\"3\" a=\"0.9\"/>";
            else if(i==p) file << "<viz:color r=\"9\" g=\"2\" b=\"200\" a=\"0.9\"/>";
            else file << "<viz:color r=\"9\" g=\"2\" b=\"2\" a=\"0.9\"/>";
            file << "</node>\n";
        }
    file << "</nodes>\n";
    file << "<edges>\n";
    for(int i = 0; i < n; ++i)
        rep(j,i)
            if(g[i][j]>=1) {
                file << "\t<edge source=\"" << i << "\" target=\"" << j << "\">";
                if(g[i][j]==2) file << "<viz:color r=\"239\" g=\"3\" b=\"3\" a=\"0.9\"/>";
                file << "</edge>\n";
            }
    file << "</edges>\n";
    file << "</graph>\n";
    file << "</gexf>\n";
}