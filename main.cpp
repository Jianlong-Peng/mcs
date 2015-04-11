/*=============================================================================
#     FileName: main.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-04-07 10:38:53
#   LastChange: 2013-04-18 19:27:39
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
#include "tools.h"
#include "mcs.h"

using namespace std;
using namespace OpenBabel;


int main(int argc, char *argv[])
{
    string prog_name(argv[0]);
    string root;
    string::size_type i = prog_name.find_last_of("/\\");
    if(i == string::npos)
	root = "./";
    else
	root = prog_name.substr(0,i+1);

    if(argc<3 || argc>9) {
	cerr << endl << "  OBJ  : to find quasi-maximum common subgraph between mol1 and mol2" << endl
	    << endl << "  Usage: " << argv[0] << "[options] mol1 mol2" << endl
            << "  [options]" << endl
            << "  --rmax n: stop recursion after 'n' steps" << endl
            << "    (default: 15000)" << endl
            << "  --smin n: eliminate those SCCSs whose cardinality < 'smin'" << endl
            << "    (default: 2)" << endl
            << "  --atom file: where atom types are defined" << endl
            << "    (default: " << root << "atom_type.txt)" << endl
            << "  --minchange: filter quasi-MCS by minimizing" << endl
            << "   \"number of mapping atoms that have different neighbors\" AND" << endl
            << "   \"number of bond changes within mapping atoms\"" << endl
            << "   (default: no such filter applied)" << endl << endl;
        exit(EXIT_FAILURE);
    }

    int rmax = 15000;
    int smin = 2;
    bool minchange = false;
    string atom_type_file = root + "atom_type.txt";
    int ii;
    for(ii=1; ii<argc; ++ii) {
        if(strncmp(argv[ii],"--",2) != 0)
            break;
        if(strcmp(argv[ii],"--rmax") == 0)
            rmax = atoi(argv[++ii]);
        else if(strcmp(argv[ii],"--smin") == 0)
            smin = atoi(argv[++ii]);
        else if(strcmp(argv[ii],"--atom") == 0)
            atom_type_file = argv[++ii];
        else if(strcmp(argv[ii],"--minchange") == 0)
            minchange = true;
        else {
            cerr << "Error: invalid option <" << argv[ii] << ">" << endl;
            exit(EXIT_FAILURE);
        }
    }
    if(argc-ii != 2) {
        cerr << "Error: invalid number of arguments, need 2 molecule files!!" << endl;
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
    
    cout << endl << "###results are based on atom types from <" << atom_type_file << ">###" << endl;
    MCSGraph *g1 = OBMol2MCSGraph(mol1,atom_type_file.c_str());
    MCSGraph *g2 = OBMol2MCSGraph(mol2,atom_type_file.c_str());
    vector<pair<vector<int>,vector<int> > > mcs_indices = findMCS(g1,g2,rmax,smin,weightOfNodeOfAG);
    for(vector<pair<vector<int>,vector<int> > >::iterator i=mcs_indices.begin(); i!=mcs_indices.end(); ++i) {
        cout << endl << "atom mapping" << endl;
        cout << " g1:";
        for(vector<int>::iterator j=i->first.begin(); j!=i->first.end(); ++j)
            printf(" %-2d",*j+1);
        cout << endl << " g2:";
        for(vector<int>::iterator j=i->second.begin(); j!=i->second.end(); ++j)
            printf(" %-2d",*j+1);
        cout << endl;
    }
    cout << endl << "totally " << mcs_indices.size() << " mappings" << endl << endl;

    if(minchange) {
        vector<int> min_indices = mapWithMinChange(g1,g2,mcs_indices);
        cout << "- atom mappings that have minimum number of mapping atoms" << endl
            << "  that have different neighbors!" << endl;
        for(vector<int>::iterator i=min_indices.begin(); i!=min_indices.end(); ++i) {
            cout << "atom mapping" << endl;
            cout << " g1:";
            for(vector<int>::iterator j=mcs_indices[*i].first.begin();
                    j!=mcs_indices[*i].first.end(); ++j)
                printf(" %-2d",*j+1);
            cout << endl << " g2:";
            for(vector<int>::iterator j=mcs_indices[*i].second.begin();
                    j!=mcs_indices[*i].second.end(); ++j)
                printf(" %-2d",*j+1);
            cout << endl << endl;
        }
    }
    delete g1;
    delete g2;
    return 0;
}

