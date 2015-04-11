/*=============================================================================
#     FileName: graph.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-04-02 16:12:17
#   LastChange: 2013-04-17 20:43:57
#      History:
=============================================================================*/

#include <iostream>
#include <utility>
#include <vector>
#include <iterator>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include "graph.h"

using std::cerr;
using std::endl;
using std::ostream;
using std::vector;
using std::pair;
using std::make_pair;

MCSMatrix::MCSMatrix(int n)
{
    assert(n > 0);
    _nnode = n;
    int _n = n*(n+1)/2;
    _m = new int [_n];
    if(_m == NULL) {
        cerr << __FILE__ << "(line " << __LINE__ << "): out of memory!" << endl;
        exit(EXIT_FAILURE);
    }
    memset(_m,0,_n*sizeof(int));
}

MCSMatrix::MCSMatrix(const MCSMatrix &m)
{
    _nnode = m.numNodes();
    int _n = _nnode*(_nnode+1)/2;
    _m = new int [_n];
    if(_m == NULL) {
        cerr << __FILE__ << "(line " << __LINE__ << "): out of memory!" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i=0; i<_n; ++i)
        for(int j=i; j<_n; ++j)
            (*this)(i,j) = m(i,j);
}

MCSMatrix::~MCSMatrix()
{
    if(_m)
        delete[] _m;
}

void MCSMatrix::expand(int n, int v)
{
    assert(n >= _nnode);
    int _n = n*(n+1)/2;
    if(_nnode == 0) {
	_m = new int [_n];
	if(_m == NULL) {
	    cerr << __FILE__ << "(line " << __LINE__ << "): out of memory!" << endl;
	    exit(EXIT_FAILURE);
	}
	_nnode = n;
	for(int i=0; i<_n; ++i)
	    _m[i] = v;
    }
    else {
	int *tmp = (int *)realloc(_m,_n);
	if(tmp == NULL) {
	    cerr << __FILE__ << "(line " << __LINE__ << "): out of memory!" << endl;
	    exit(EXIT_FAILURE);
	}
	_m = tmp;
	for(int i=_nnode*(_nnode+1)/2; i<_n; ++i)
	    _m[i] = v;
        _nnode = n;
    }
}

MCSMatrix &MCSMatrix::operator=(const MCSMatrix &m)
{
    if(this == &m)
        return *this;
    if(_nnode==0 && m.numNodes()!=0) {
	_m = new int [m.numNodes()];
	if(_m == NULL) {
	    cerr << __FILE__ << "(line " << __LINE__ << "): out of memory!" << endl;
	    exit(EXIT_FAILURE);
	}
	_nnode = m.numNodes();
    }
    assert(_nnode==m.numNodes());
    for(int i=0; i<_nnode; ++i)
        for(int j=i; j<_nnode; ++j)
	    (*this)(i,j) = m(i,j);
    return *this;
}

int MCSMatrix::getIndex(int i, int j) const
{
    assert(_nnode!=0 && i>=0 && j>=0 && i<_nnode && j<_nnode);
    int index;
    if(i < j)
        index = (2*_nnode-i+1)*i/2 + j-i;
    else
        index = (2*_nnode-j+1)*j/2 + i-j;
    return index;
}

const int &MCSMatrix::operator()(int i, int j) const
{
    return _m[getIndex(i,j)];
}

int &MCSMatrix::operator()(int i, int j)
{
    return _m[getIndex(i,j)];
}

void MCSMatrix::display(ostream &fout, const char *sep) const
{
    if(_nnode == 0)
	fout << "empty matrix" << endl;
    else {
        for(int i=0; i<_nnode; ++i) {
	    for(int j=0; j<_nnode; ++j)
		fout << sep << (*this)(i,j);
	    fout << endl;
	}
    }
}

bool MCSGraph::connect(int i, int j) const 
{
    if(i==j || _m(i,j)==0)
        return false;
    else
        return true;
}

const int &MCSGraph::getNodeLabelOf(int i) const
{
    assert(i>=0 && i<static_cast<int>(_vertex.size()));
    return _vertex[i];
}

int &MCSGraph::getNodeLabelOf(int i)
{
    assert(i>=0 && i<static_cast<int>(_vertex.size()));
    return _vertex[i];
}

void MCSGraph::display(ostream &fout, const char *sep) const
{
    fout << "vertex labels:";
    for(vector<int>::const_iterator i=_vertex.begin(); i!=_vertex.end(); ++i)
        fout << sep << *i;
    fout << endl;
    _m.display(fout,sep);
}

vector<int> MCSGraph::getNborsOf(int i) const
{
    int num_nodes = _m.numNodes();
    assert(i>=0 && i<num_nodes);
    vector<int> nbors;
    for(int j=0; j<num_nodes; ++j) {
        if(i==j || _m(i,j)==0)
            continue;
        nbors.push_back(j);
    }
    return nbors;
}

static float nodeValue(const MCSGraph *g1, int i1, const MCSGraph *g2, int i2)
{
    if(g1->getNodeLabelOf(i1) == g2->getNodeLabelOf(i2))
        return 1.;
    else
        return 0.;
}
AssocGraph::AssocGraph(const MCSGraph *g1, const MCSGraph *g2,
            float (*nodeOfAG)(const MCSGraph *g1, int i1, const MCSGraph *g2, int i2),
            int edge_type)
{
    assert(edge_type==1 || edge_type==2);
    int node1 = g1->numNodes();
    int node2 = g2->numNodes();
    // figure out the number of nodes to be created!
    // fill '_node' and '_indices'
    int num = 0;
    for(int i=0; i<node1; ++i) {
	for(int j=0; j<node2; ++j) {
	    float value;
	    if(nodeOfAG)
		value = nodeOfAG(g1,i,g2,j);
	    else
		value = nodeValue(g1,i,g2,j);
	    if(value == 0.)
		continue;
	    ++num;
	    _node.push_back(value);
	    _indices.push_back(make_pair(i,j));
	}
    }
    _m.expand(num);
    // fill edge of AG
    for(int i=0; i<num-1; ++i) {
	int index_1i = _indices[i].first;
	int index_2i = _indices[i].second;
	for(int j=i+1; j<num; ++j) {
	    int index_1j = _indices[j].first;
	    int index_2j = _indices[j].second;
	    if(index_1i==index_1j || index_2i==index_2j)
		continue;
	    bool connect_1i_1j = g1->connect(index_1i,index_1j);
	    bool connect_2i_2j = g2->connect(index_2i,index_2j);
            if(edge_type == 1) {
                if(connect_1i_1j == connect_2i_2j)
                    _m(i,j) = 1;
            }  else if(edge_type == 2) {
                if(connect_1i_1j==true && connect_2i_2j==true)
                    _m(i,j) = 1;
            } else
                ;
	}
    }
}

pair<vector<int>,vector<int> > AssocGraph::toGraphIndices(const vector<int> &index) const
{
    vector<int> index1;
    vector<int> index2;
    for(vector<int>::const_iterator i=index.begin(); i!=index.end(); ++i) {
	const pair<int,int> &temp = toGraphIndex(*i);
	index1.push_back(temp.first);
	index2.push_back(temp.second);
    }
    return make_pair(index1,index2);
}

int AssocGraph::toAGIndex(int i, int j) const
{
    for(vector<pair<int,int> >::size_type k=0; k<_indices.size(); ++k) {
        if(_indices[k].first==i && _indices[k].second==j)
            return static_cast<int>(k);
    }
    return -1;
}

vector<int> AssocGraph::toAGIndex(const vector<int> &index_i, const vector<int> &index_j) const
{
    if(index_i.size() != index_j.size()) {
        cerr << __FILE__ << "(line " << __LINE__ << ") warnning: 1st and 2nd parameter have different number of indices!!!" << endl;
        return vector<int>();
    }
    vector<int> result;
    for(vector<int>::size_type k=0; k<index_i.size(); ++k)
        result.push_back(toAGIndex(index_i[k],index_j[k]));
    return result;
}

vector<int> AssocGraph::getNborsOf(int i) const
{
    int num_nodes = _m.numNodes();
    assert(i>=0 && i<num_nodes);
    vector<int> nbors;
    for(int j=0; j<num_nodes; ++j) {
        if(i==j || _m(i,j)==0)
            continue;
        nbors.push_back(j);
    }
    return nbors;
}


