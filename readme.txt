
Usage
=====
// 1. prepare two graphs g1 and g2
MCSGraph *g1 = ...
MCSGraph *g2 = ...
// 2. construct association graph from g1 and g2
AssocGraph ag(g1,g2,nodeOfAG);
// 3. find all cliques
vector<vector<int> > cliques = findCliques(ag.getMatrix());
// 4. get the maximum clique (MCL)
vector<int> mcl = maxClique(cliques,ag.getNodeValues());
// 5. get atom mapping of the MCL
pair<vector<int>,vector<int> > graph_index = ag.toGraphIndex(mcl);
