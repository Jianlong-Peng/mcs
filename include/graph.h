/*=============================================================================
#     FileName: graph.h
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-04-02 16:07:56
#   LastChange: 2013-05-09 04:06:06
#      History:
=============================================================================*/

#ifndef  GRAPH_H
#define  GRAPH_H

#include <iostream>
#include <vector>
#include <utility>
#include <cassert>

// upper triangle square matrix with diagonal elements kept
class MCSMatrix
{
public:
    MCSMatrix():_m(0),_nnode(0) {}
    explicit MCSMatrix(int n);
    MCSMatrix(const MCSMatrix &m);
    ~MCSMatrix();
    // to expand the matrix. 'n' must be >= '_nnode'
    // newly elements will be set to 'v'
    void expand(int n, int v=0);
    // if this->_nnode==0, allocate memory first!!
    // otherwise, this->_nnode must be equal to m.numNodes()
    MCSMatrix &operator=(const MCSMatrix &m);
    const int &operator()(int i, int j) const;
    int &operator()(int i, int j);
    int numNodes() const {return _nnode;}
    void display(std::ostream &fout, const char *sep) const;
private:
    int getIndex(int i, int j) const;
    int *_m; // of length n*(n+1)/2
    int _nnode;
};

class MCSGraph
{
public:
    explicit MCSGraph(const std::vector<int> &vertex)
        :_vertex(vertex),_m(static_cast<int>(_vertex.size()))
    {}
    MCSGraph(const std::vector<int> &vertex, const MCSMatrix &m)
        :_vertex(vertex),_m(m)
    {}
    MCSGraph(const int vertex[], int n): _vertex(vertex,vertex+n),_m(n) {}
    ~MCSGraph() {};
    // get edge
    const int &operator()(int i, int j) const {return _m(i,j);}
    int &operator()(int i, int j) {return _m(i,j);}
    // if two nodes are connected
    bool connect(int i, int j) const; // false if (i,j)==0 or i==j
    // number of nodes
    int numNodes() const {return _m.numNodes();}
    // get node labels
    const std::vector<int> &getNodeLabels() const {return _vertex;}
    std::vector<int> &getNodeLabels() {return _vertex;}
    // get node i's label
    const int &getNodeLabelOf(int i) const;
    int &getNodeLabelOf(int i);
    // display nodes and edges
    void display(std::ostream &fout, const char *sep) const;
    // get neighbors
    std::vector<int> getNborsOf(int i) const;
private:
    std::vector<int> _vertex; // vector of vertexes
    MCSMatrix _m;
};

// for association graph
class AssocGraph
{
public:
    // Parameter:
    // - nodeOfAG: to determine the node value(weight) of association graph
    //     default:
    //       1.0 - g1->getNodeLabelOf(i1)==g2->getNodeLabelOf(i2)
    //       0.0 - otherwise
    //     if 'nodeOfAG' returns 0, the corresponding node will not be created.
    // - edge_type: the way to make edge <default: 1>
    //     1 - e(1i2i,1j2j) equals to 1 if it satisfies one of following conditions:
    //         1) connect(1i,1j)==true and connect(2i,2j)==true
    //         2) connect(1i,1j)==false and connect(2i,2j)==false
    //     2 - e(1i2i,1j2j) equals to 1 if it satisfies that
    //         "connect(1i,1j)==true and connect(2i,2j)==true"
    AssocGraph(const MCSGraph *g1, const MCSGraph *g2,
            float (*nodeOfAG)(const MCSGraph *g1, int i1, const MCSGraph *g2, int i2)=NULL,
            int edge_type=1);
    const std::vector<std::pair<int,int> > &getPairIndices() const {return _indices;}
    std::vector<std::pair<int,int> > &getPairIndices() {return _indices;}
    // get corresponding indices of 'g1' and 'g2'
    const std::pair<int,int> &toGraphIndex(int i) const {assert(i<_node.size());return _indices[i];}
    std::pair<int,int> &toGraphIndex(int i) {assert(i<_node.size());return _indices[i];}
    std::pair<std::vector<int>,std::vector<int> > toGraphIndices(const std::vector<int> &index) const;
    // graph pair index ==> index of AG
    int toAGIndex(int i, int j) const;
    std::vector<int> toAGIndex(const std::vector<int> &index_i, const std::vector<int> &index_j) const;
    // matrix
    const MCSMatrix &getMatrix() const {return _m;}
    MCSMatrix &getMatrix() {return _m;}
    // node values
    const std::vector<float> getNodeValues() const {return _node;}
    std::vector<float> getNodeValues() {return _node;}
    float getNodeValueOf(int i) const {return _node[i];}
    // return nbors
    std::vector<int> getNborsOf(int i) const;
private:
    std::vector<float> _node; // node values
    std::vector<std::pair<int,int> > _indices; // original graphs' indices
    MCSMatrix _m;
};

/*
// OBJ: to create a association graph from two graphs.
// connect: used to determine the edge(1 or 0) of association graph
//          if NULL, using member function 'connect' of MCSGraph
// e(1i2i,1j2j)=1 if only if it follows one of two conditions:
//  1) connect(1i,1j)=true and connect(2i,2j)=true
//  2) connect(1i,2j)=false and connect(2i,2j)=false
// However, if 1i==1j or 2i==2j, here we let e(1i2i,1j2j) be 0!
MCSMatrix *makeAG(const MCSGraph &g1, const MCSGraph &g2,
        bool (*connect)(const MCSGraph &g, int i, int j)=NULL);

// OBJ: convert index of an association graph to indices of graphs.
//      suppose that the association graph was constructed from 'makeAG'.
// Parameter
//   - nnodes2: number of nodes of graph 'g2'
//   - agindex: index of a association graph
//   - i,j: index of graph 1 and 2, respectively
void AGIndex2GraphIndex(int nnodes2, int agindex, int *i, int *j);
*/
#endif   /* ----- #ifndef GRAPH_H  ----- */

