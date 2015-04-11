/*=============================================================================
#     FileName: atom_type.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-04-08 19:53:17
#   LastChange: 2013-04-08 20:26:52
#      History:
=============================================================================*/
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <openbabel/atom.h>
#include "atom_type.h"

using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;

static vector<string> _atom_types;

static vector<string> split(const string &s)
{
    vector<string> ss;
    string::size_type i=0;
    while(i < s.size()) {
        while(i<s.size() && (s[i]==' ' || s[i]=='\t')) ++i;
        if(i == s.size()) break;
        string::size_type j = s.find_first_of(" \t",i);
        if(j != string::npos) {
            ss.push_back(s.substr(i,j-i));
            i = j+1;
        } else {
            ss.push_back(s.substr(i));
            i = s.size();
        }
    }
    return ss;
}

static void read_atom_types(const char *atomtype_file)
{
    ifstream inf(atomtype_file);
    if(!inf) {
        cerr << __FILE__ << "(line " << __LINE__ << "): failed to open file " << atomtype_file << endl;
        exit(EXIT_FAILURE);
    }
    string line;
    while(getline(inf,line)) {
        if(line.empty() || line[0]=='#')
            continue;
        vector<string> temp = split(line);
        _atom_types.push_back(temp[0]);
    }
    inf.close();
}

int atom_type(OpenBabel::OBAtom &atom, const char *atomtype_file)
{
    if(_atom_types.empty())
        read_atom_types(atomtype_file);
    vector<string>::size_type type;
    for(vector<string>::size_type i=0; i<_atom_types.size(); ++i) {
        if(atom.MatchesSMARTS(_atom_types[i].c_str()))
            type = i+1;
    }
    return static_cast<int>(type);
}

