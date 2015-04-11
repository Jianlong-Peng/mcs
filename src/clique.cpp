/*=============================================================================
#     FileName: clique.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-04-03 17:12:55
#   LastChange: 2013-04-10 19:10:58
#      History:
=============================================================================*/

#include <vector>
#include <iostream>
#include <utility>
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include "graph.h"
#include "clique.h"

using std::vector;
using std::cerr;
using std::endl;
using std::make_pair;
using std::max_element;

#if !defined(EXTEND_VERSION1)
#define EXTEND_VERSION1
#endif

static int Rmax = 0;
static int Steps = 0;

static void extend1(const MCSMatrix &m, int *old, int ne, int ce,
        vector<int> &compsub, vector<vector<int> > &cliques/*, int times*/)
{
    ++Steps;
    bool sucexp;
    do {
	/*
	std::cout << endl << "ne=" << ne << ", old: ";
	copy(old,old+ce+1,std::ostream_iterator<int>(std::cout," "));
	std::cout << endl;
	*/
        // 1. to see if there exists node in old[0:ne] such that it connects to
        //    each node in old[ne+1:ce]
        sucexp = true;
        for(int i=0; i<=ne && sucexp; ++i) {
            bool allcon = true;
            for(int j=ne+1; j<=ce && allcon; ++j) {
                if(m(old[i],old[j]) == 0)
                    allcon = false;
                else
                    allcon = true;
            }
            sucexp = !allcon;
        }
        // 2. select candidate
        if(sucexp) {
            int sel = old[ne+1];
            int newne = -1;
            int *new_ = new int [ce];
            // 2.1 fill new set 'not' by removing all nodes, which are not connected to 'sel', 
            //     in old[0:ne]
            int i = 0;
            for(; i<=ne; ++i) {
                if(m(sel,old[i]))
                    new_[++newne] = old[i];
            }
            // 2.2 fill new set 'candidate' by removing all node, which are not connected to sel,
            //     in old[ne+1:ce]
            int newce = newne;
            for(i+=1; i<=ce; ++i) {
                if(m(sel,old[i]))
                    new_[++newce] = old[i];
            }
	    /*
	    std::cout << "sel=" << sel << ", ne=" << ne << ", newne=" << newne << endl;
	    std::cout << "new_: ";
	    copy(new_,new_+newce+1,std::ostream_iterator<int>(std::cout," "));
	    std::cout << endl;
	    */
            // 2.3 add 'sel' to compsub
            compsub.push_back(sel);
	    /*
	    std::cout << "compsub: ";
	    copy(compsub.begin(),compsub.end(),std::ostream_iterator<int>(std::cout," "));
	    std::cout << endl;
	    */
            // 2.3.1. find a clique
            if(newce==-1 || (Rmax>0 && Steps>=Rmax)) {
                cliques.push_back(compsub);
		//std::cout << "find a clique" << endl;
	    }
            // 2.3.2. go on extending the current complete subgraph
            else if(newne < newce)
                extend1(m,new_,newne,newce,compsub,cliques/*,times+1*/);
            // 2.3.3. no more candidates, but the current 'sel' have been involved 
            //        in previous complete subgraph
            else
                ;
            compsub.pop_back();
            ++ne;
            delete[] new_;
        }
    }while(sucexp);
}

static void extend2(const MCSMatrix &m, int *old, int ne, int ce,
        vector<int> &compsub, vector<vector<int> > &cliques/*, int times*/)
{
    int minnod = ce+1;  // minnod: minimum number of disconnect
    int nod = 0, fixp=-1, s=-1;
    ++Steps;
    // 1. to determine number of disconnect nodes in old[ne+1:ce] for each node in 'old',
    //    and look for the node with minimum 'minnod'.
    for(int i=0; i<=ce && minnod!=0; ++i) {
        int pos=-1;
        int count = 0;
        // 1.1 count disconnections
        for(int j=ne+1; j<=ce; ++j) {
            if(i == j) // in case that (*m)(old[i],old[i]) = 0 !!!!!!!!!
                continue;
            if(m(old[i],old[j]) == 0) {
                ++count;
                pos = j;
            }
        }
        // 1.2 test new minimum
        if(count < minnod) {
            fixp = old[i];  // the node with minimum number of disconnect nodes in old[ne+1:ce]
            minnod = count;
            if(i <= ne) s=pos; // node in old[ne+1:ce] which was disconnect to 'fixp'
            else {
                s = i;
                nod = 1;
                // whenever 'fixp' was in old[ne+1:ce], it'll first check this node and
                // move the 'fixp' to 'not'. So the while-loop in step 2 needs to be run
                // one more time.
            }
        }
    }
    nod += minnod;
    // 2. select candidates
    int *new_ = new int [ce+1];
    while(nod > 0) {
        // 2.1 interchange 'sel' and old[ne+1]
        int sel = old[s];
        old[s] = old[ne+1];
        old[ne+1] = sel;
        // 2.2 fill new set 'not' by removing nodes disconnected to 'sel' in old[0:ne]
        int newne=-1;
        /*
        int *new_ = new int [ce+1];
        if(new_ == NULL) {
            cerr << __FILE__ << "(line " << __LINE__ << "): out of memory!" << endl;
            exit(EXIT_FAILURE);
        }
        */
        int i=0;
        for(; i<=ne; ++i) {
            if(m(sel,old[i]))
                new_[++newne] = old[i];
        }
        // 2.3 fill new set 'candidate' by removing nodes disconnected to 'sel' in old[ne:ce-1]
        int newce = newne;
        for(i+=1; i<=ce; ++i) {
            if(m(sel,old[i]))
                new_[++newce] = old[i];
        }
        // add 'sel' to 'compsub'
        compsub.push_back(sel);
        if(newce==-1 || (Rmax>0 && Steps>=Rmax)) {
            cliques.push_back(compsub);
	    /*
	    std::cout << endl << "compsub: ";
	    copy(compsub.begin(),compsub.end(),std::ostream_iterator<int>(std::cout," "));
	    std::cout << endl;
	    */
	}
        else if(newne < newce)
            extend2(m,new_,newne,newce,compsub,cliques/*,times+1*/);
        else
            ;
        compsub.pop_back();
        ++ne;
        --nod;
        //delete[] new_;
        // look for another node that is disconnect to 'fixp'
        if(nod > 0) {
            for(s=ne+1; m(fixp,old[s]); ++s)
                ;
        }
    }
    delete[] new_;
}

vector<vector<int> > findCliques(const MCSMatrix &m, int rmax)
{
    int nnodes = m.numNodes();
    int *old = new int [nnodes];
    if(old == NULL) {
        cerr << __FILE__ << "(line " << __LINE__ << "): out of memory!" << endl;
        exit(EXIT_FAILURE);
    }
    for(int i=0; i<nnodes; ++i)
        old[i] = i;
    vector<vector<int> > cliques;
    vector<int> compsub;
    Rmax = rmax;
    Steps = 0;
#if defined(EXTEND_VERSION1)
    extend1(m,old,-1,nnodes-1,compsub,cliques);
#else
    extend2(m,old,-1,nnodes-1,compsub,cliques);
#endif
    cerr << "Totally " << Steps << " recursion steps!" << endl;
    delete[] old;
    return cliques;
}

vector<vector<int> > maxClique(const vector<vector<int> > &cliques,
        const vector<float> &node_weight)
{
    vector<float> sums;
    vector<vector<int> > result;
    for(vector<vector<int> >::size_type i=0; i<cliques.size(); ++i) {
        float sum = 0.;
        for(vector<int>::const_iterator j=cliques[i].begin(); j!=cliques[i].end(); ++j)
            sum += node_weight[*j];
	sums.push_back(sum);
    }
    float max_sum = *max_element(sums.begin(),sums.end());
    for(vector<vector<int> >::size_type i=0; i<cliques.size(); ++i) {
	if(sums[i] == max_sum)
	    result.push_back(vector<int>(cliques[i].begin(),cliques[i].end()));
    }
    return result;
}


