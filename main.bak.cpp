/*=============================================================================
#     FileName: main.bak.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-04-07 10:38:53
#   LastChange: 2013-04-18 19:27:23
#      History:
=============================================================================*/
#include <iostream>
#include <vector>
#include <utility>
#include <iterator>
#include <string>
#include <cstdio>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include "graph.h"
#include "clique.h"
#include "tools.h"

using namespace std;
using namespace OpenBabel;

vector<pair<vector<int>,vector<int> > > mcs(OBMol &mol1, OBMol &mol2, int rmax);

static string atom_type_file;

int main(int argc, char *argv[])
{
    string prog_name(argv[0]);
    string root;
    string::size_type i = prog_name.find_last_of("/\\");
    if(i == string::npos)
	root = "./";
    else
	root = prog_name.substr(0,i+1);

    if(argc<3 || argc>7) {
	cerr << endl << "  OBJ  : to find maximum common subgraph between mol1 and mol2" << endl
	    << endl << "  Usage: " << argv[0] << "[options] mol1 mol2" << endl
            << "  [options]" << endl
            << "  --rmax n: stop recursion after 'n' steps" << endl
            << "    (default: do not stop)" << endl
            << "  --atom file: where atom types are defined" << endl
            << "    (default: " << root << "atom_type.txt)" << endl << endl
            << "  Attention: for molecules with lots of atoms, it could be very time-consuming!" << endl
            << "             In that case, you can use './mcs' to get a quasi-MCS instead!" << endl << endl;
        exit(EXIT_FAILURE);
    }

    int rmax = 0;
    atom_type_file = root + "atom_type.txt";
    int ii;
    for(ii=1; ii<argc; ) {
        if(strncmp(argv[ii],"--",2) != 0)
            break;
        if(strcmp(argv[ii],"--rmax") == 0)
            rmax = atoi(argv[ii+1]);
        else if(strcmp(argv[ii],"--atom") == 0)
            atom_type_file = argv[ii+1];
        else {
            cerr << "Error: invalid option <" << argv[ii] << ">" << endl;
            exit(EXIT_FAILURE);
        }
        ii += 2;
    }
    if(argc-ii != 2) {
        cerr << "Error: invalid number of arguments, need 2 molecule files" << endl;
        exit(EXIT_FAILURE);
    }
    string name1(argv[ii]);
    string name2(argv[ii+1]);
    string format1 = name1.substr(name1.rfind(".")+1);
    string format2 = name2.substr(name2.rfind(".")+1);
    OBConversion conv;
    OBMol mol1,mol2;
    if(!conv.SetInFormat(format1.c_str()) || !conv.ReadFile(&mol1,name1)) {
        cerr << "Error: failed to read molecule from " << name1 << endl;
        exit(EXIT_FAILURE);
    }
    if(!conv.SetInFormat(format2.c_str()) || !conv.ReadFile(&mol2,name2)) {
        cerr << "Error: failed to read molecule from " << name2 << endl;
        exit(EXIT_FAILURE);
    }
    vector<pair<vector<int>,vector<int> > > graph_indices = mcs(mol1,mol2,rmax);

    for(vector<pair<vector<int>,vector<int> > >::iterator i=graph_indices.begin(); i!=graph_indices.end(); ++i) {
        cout << endl << "atom mapping:" << endl;
	cout << " g1:";
        for(vector<int>::iterator j=(*i).first.begin();j!=(*i).first.end();++j)
	    printf(" %-2d",*j+1);
        cout << endl << " g2:";
	for(vector<int>::iterator j=(*i).second.begin();j!=(*i).second.end();++j)
	    printf(" %-2d",*j+1);
        cout << endl;
    }
    cout << endl << "totally " << graph_indices.size() << " atom mappings for MCS" << endl << endl;

    return 0;
}

vector<pair<vector<int>,vector<int> > > mcs(OBMol &mol1, OBMol &mol2, int rmax)
{
    MCSGraph *g1 = OBMol2MCSGraph(mol1,atom_type_file.c_str());
    MCSGraph *g2 = OBMol2MCSGraph(mol2,atom_type_file.c_str());
    //cerr << "to make association graph" << endl;
    AssocGraph ag(g1,g2,weightOfNodeOfAG);
    //cerr << "to find cliques" << endl;
    vector<vector<int> > cliques = findCliques(ag.getMatrix(),rmax);
    cerr << "-> there are " << cliques.size() << " cliques" << endl;
    /*
    for(vector<vector<int> >::iterator i=cliques.begin(); i!=cliques.end(); ++i) {
	if(i->size() > 10) {
	    float sum = 0;
	    for(vector<int>::iterator j=i->begin(); j!=i->end(); ++j)
		sum += ag.getNodeValueOf(*j);
	    cout << "clique: ";
	    copy(i->begin(),i->end(),std::ostream_iterator<int>(cout," "));
	    cout << endl;
	    cout << "atom mapping with " << i->size() << " maps; weight=" << sum << endl;
	    pair<vector<int>,vector<int> > temp = ag.toGraphIndex(*i);
	    for(vector<int>::iterator j=temp.first.begin(); j!=temp.first.end(); ++j)
		cout << " " << *j+1;
	    cout << endl;
	    for(vector<int>::iterator j=temp.second.begin(); j!=temp.second.end(); ++j)
		cout << " " << *j+1;
	    cout << endl;
	}
    }
    */
    vector<vector<int> > mcl = maxClique(cliques,ag.getNodeValues());
    cerr << "-> there are " << mcl.size() << " maximum cliques" << endl;
    //cout << "MCL: ";
    //copy(mcl.begin(),mcl.end(),ostream_iterator<int>(cout," "));
    //cerr << "to convert clique index to graph index" << endl;
    vector<pair<vector<int>,vector<int> > > graph_indices;
    for(vector<vector<int> >::iterator i=mcl.begin(); i!=mcl.end(); ++i)
	graph_indices.push_back(ag.toGraphIndices(*i));
    delete g1;
    delete g2;
    return graph_indices;
}

