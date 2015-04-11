/*=============================================================================
#     FileName: tools.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-04-08 19:53:17
#   LastChange: 2013-04-11 10:41:32
#      History:
=============================================================================*/
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <cstdlib>
#include <openbabel/atom.h>
#include <openbabel/mol.h>
#include "tools.h"

using std::cerr;
using std::endl;
using std::vector;
using std::map;
using std::string;
using std::ifstream;
using namespace OpenBabel;

static map<string,int> _atom_types;
static const int STEP = 1000;
static int MAX_TYPE = 0; // for default atom type

static void read_atom_types(const char *atomtype_file)
{
    ifstream inf(atomtype_file);
    if(!inf) {
        cerr << __FILE__ << "(line " << __LINE__ << "): failed to open file " << atomtype_file << endl;
        exit(EXIT_FAILURE);
    }
    string line;
    int prefix = -STEP;
    int i = -1;
    while(getline(inf,line)) {
	if(i == -1) { // begining block
	    prefix += STEP;
	    i = 0;
	}
	else if(line == "") { // new block
	    prefix += STEP;
	    i = 0;
	    getline(inf,line);
	}
	while(!inf.eof() && (line=="" || line[0]=='#'))
	    getline(inf,line);
	if(inf.eof())
	    break;
	string::size_type j = line.find_first_of(" \t");
	if(j==string::npos)
	    _atom_types[line] = prefix + i;
	else
	    _atom_types[line.substr(0,j)] = prefix + i;
	++i;
    }
    if(prefix==-STEP || (prefix==0 && i==0)) {
	cerr << atomtype_file << ": contains no pattern!" << endl;
	exit(EXIT_FAILURE);
    }
    MAX_TYPE = prefix + STEP;
}    

int atomType(OBAtom &atom, const char *atomtype_file)
{
    if(_atom_types.empty())
        read_atom_types(atomtype_file);
    int type = MAX_TYPE;
    for(map<string,int>::iterator i=_atom_types.begin(); i!=_atom_types.end(); ++i) {
        if(atom.MatchesSMARTS((i->first).c_str()))
            type = i->second;
    }
    return type;
}

float weightOfNodeOfAG(const MCSGraph *g1, int i1, const MCSGraph *g2, int i2)
{
    int label1 = g1->getNodeLabelOf(i1);
    int label2 = g2->getNodeLabelOf(i2);
    if(label1 == label2)
	return 1.;
    else {
	int prefix1 = label1 / STEP;
	int prefix2 = label2 / STEP;
	if(prefix1 == prefix2)
	    return 0.5;
	else
	    return 0.;
    }
}

MCSGraph *OBMol2MCSGraph(OBMol &mol, const char *atom_type_file)
{
    vector<int> labels;
    FOR_ATOMS_OF_MOL(atom,mol)
	labels.push_back(atomType(*atom,atom_type_file));
#ifdef DEBUG
    cout << "atom types" << endl << " ";
    copy(labels.begin(),labels.end(),std::ostream_iterator<int>(cout," "));
    cout << endl;
#endif
    MCSGraph *g = new MCSGraph(labels);
    FOR_BONDS_OF_MOL(bond,mol) {
        int i = bond->GetBeginAtom()->GetIdx();
        int j = bond->GetEndAtom()->GetIdx();
        int type = 6;
        if(bond->IsSingle())
            type = 1;
        if(bond->IsDouble())
            type = 2;
        if(bond->IsTriple())
            type = 3;
        if(bond->IsAromatic())
            type = 4;
        if(bond->IsAmide())
            type = 5;
        (*g)(i-1,j-1) = type;
    }
    return g;
}

