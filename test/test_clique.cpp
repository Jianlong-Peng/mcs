#include <iostream>
#include <vector>
#include <iterator>
#include "graph.h"
#include "clique.h"

using namespace std;

int main()
{
    int labels[] = {1,1,1,1,1,1};
    MCSGraph g1(labels,3);
    MCSGraph g2(labels,4);
    
    // add edge
    // g1: 0-1-2
    // g2: 0-1-3
    //       |
    //       2
    g1(0,1) = 1;
    g1(1,2) = 1;
    g2(0,1) = 1;
    g2(1,2) = 1;
    g2(1,3) = 1;
    // display g1 and g2
    cout << "g1" << endl;
    g1.display(cout," ");
    cout << "g2" << endl;
    g2.display(cout," ");
    
    // association graph
    AssocGraph ag(&g1,&g2);
    cout << "\nassociation graph" << endl;
    ag.getMatrix().display(cout," ");
    cout << endl;
    
    // find cliques
    vector<vector<int> > cliques = findCliques(ag.getMatrix());
    for(vector<vector<int> >::size_type i=0; i<cliques.size(); ++i) {
        cout << "clique indices: ";
        copy(cliques[i].begin(),cliques[i].end(),ostream_iterator<int>(cout," "));
        cout << endl;
        pair<vector<int>,vector<int> > graph_index = ag.toGraphIndex(cliques[i]);
        cout << "graph mapping" << endl;
        cout << "  ";
        copy(graph_index.first.begin(),graph_index.first.end(),ostream_iterator<int>(cout," "));
        cout << endl << "  ";
        copy(graph_index.second.begin(),graph_index.second.end(),ostream_iterator<int>(cout," "));
        cout << endl << endl;
    }
    cout << "totally " << cliques.size() << " cliques" << endl << endl;
    
    // maximum clique
    vector<vector<int> > mcl = maxClique(cliques,ag.getNodeValues());
    for(vector<vector<int> >::iterator i=mcl.begin(); i!=mcl.end(); ++i) {
        pair<vector<int>,vector<int> > graph_index = ag.toGraphIndex(*i);
        cout << "maximum clique: ";
	copy((*i).begin(),(*i).end(),ostream_iterator<int>(cout," "));
        cout << endl << "graph mapping" << endl;
        cout << "  ";
        copy(graph_index.first.begin(),graph_index.first.end(),ostream_iterator<int>(cout," "));
        cout << endl << "  ";
        copy(graph_index.second.begin(),graph_index.second.end(),ostream_iterator<int>(cout," "));
	cout << endl << endl;
    }
    
    return 0;
}
