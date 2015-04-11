/*=============================================================================
#     FileName: mcs.cpp
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-04-09 16:37:41
#   LastChange: 2013-04-18 19:31:30
#      History:
=============================================================================*/
#include <iostream>
#include <vector>
#include <deque>
#include <bitset>
#include <utility>
#include <queue>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <cassert>
#include "graph.h"
#include "clique.h"

using std::vector;
using std::deque;
using std::bitset;
using std::pair;
using std::queue;
using std::cout;
using std::cerr;
using std::endl;

// used for bitset
const int NUMNBORS = 10;
//#define DEBUG

static void print_match(const pair<vector<int>,vector<int> > &match)
{
    cout << "  g1: ";
    for(vector<int>::const_iterator jj=match.first.begin(); jj!=match.first.end(); ++jj)
        cout << " " << *jj+1;
    cout << endl << "  g2: ";
    for(vector<int>::const_iterator jj=match.second.begin(); jj!=match.second.end(); ++jj)
        cout << " " << *jj+1;
    cout << endl;
}

// find all connect graphs; called by 'findMCS'
// return: vector(vector(index of nodes))
static vector<vector<int> > connectGraph(const MCSGraph *g, const vector<int> &nodes)
{
    vector<vector<int> > result;
    int num_nodes = static_cast<int>(nodes.size());
    vector<bool> visited(num_nodes,false);
    int j=0;
    queue<int> toBeVisit;
    while(1) {
        int i = j;
        while(i<num_nodes && visited[i])
            ++i;
        if(i == num_nodes)
            break;
        j = i+1;
        toBeVisit.push(i);
        visited[i] = true;
        vector<int> graph;
        while(!toBeVisit.empty()) {
            int k = toBeVisit.front();
            graph.push_back(k);
            for(int m=0; m<num_nodes; ++m) {
                if(m==k || visited[m])
                    continue;
                if((*g)(nodes[k],nodes[m])) {
                    toBeVisit.push(m);
                    visited[m] = true;
                }
            }
            toBeVisit.pop();
        }
        result.push_back(graph);
    }
    return result;
}
// static method called by 'extend'
// to see if any of graph id of 'nbor_id' has already in group
static bool compatible(const AssocGraph &ag, const vector<int> &group, int nbor_id)
{
    pair<vector<int>,vector<int> > pair_indices = ag.toGraphIndices(group);
    pair<int,int> pair_index = ag.toGraphIndex(nbor_id);
    if(find(pair_indices.first.begin(),pair_indices.first.end(),pair_index.first) !=
            pair_indices.first.end())
        return false;
    if(find(pair_indices.second.begin(),pair_indices.second.end(),pair_index.second) !=
            pair_indices.second.end())
        return false;
    return true;
}
// static method called by 'extend'
static void extendGroup(const AssocGraph &ag, vector<vector<int> > &groups, vector<int> &candidate)
{
    if(candidate.empty())
        return;
    if(groups.empty()) {
        for(vector<int>::iterator i=candidate.begin(); i!=candidate.end(); ++i)
            groups.push_back(vector<int>(1,*i));
        return;
    }
    vector<vector<int> > new_groups;
    for(vector<vector<int> >::iterator i=groups.begin(); i!=groups.end(); ++i) {
        for(vector<int>::iterator j=candidate.begin(); j!=candidate.end(); ++j) {
            vector<int> temp(i->begin(),i->end());
            if(compatible(ag,temp,*j)) {
                temp.push_back(*j);
                new_groups.push_back(temp);
            }
        }
    }
    if(new_groups.empty())
        return ;
    groups.clear();
    for(vector<vector<int> >::iterator i=new_groups.begin(); i!=new_groups.end(); ++i)
        groups.push_back(*i);
}
// static method(recursion) called by 'extendSCCS'
static void extend(const AssocGraph &ag, const MCSGraph *g1, const MCSGraph *g2,
        vector<bool> &visited1, vector<bool> &visited2,
        vector<vector<int> > &result, vector<int> &temp,
        deque<int> &nborsToBeGet)
{
#ifdef DEBUG
    cout << " in 'extend', nborsToBeGet: ";
    for(deque<int>::iterator ii=nborsToBeGet.begin(); ii!=nborsToBeGet.end(); ++ii)
        cout << " " << *ii + 1;
    cout << endl;
#endif
    vector<int> deleted_ids;
    int curr_id;
    pair<int,int> curr_pair_index;
    vector<int> nbors;
    do {
        curr_id = nborsToBeGet.front();
        nborsToBeGet.pop_front();
        deleted_ids.push_back(curr_id);
        curr_pair_index = ag.toGraphIndex(curr_id);
#ifdef DEBUG
        cout << " to check neighbor matching of " << curr_id+1 << " (" << curr_pair_index.first+1 << "," << curr_pair_index.second+1 << ") " << endl;
#endif
        // remove invalid neighbors
        nbors.clear();
        nbors = ag.getNborsOf(curr_id);
        for(vector<int>::iterator i=nbors.begin(); i!=nbors.end();) {
            pair<int,int> nbor_pair_index = ag.toGraphIndex(*i);
            if(visited1[nbor_pair_index.first] || visited2[nbor_pair_index.second])
                i = nbors.erase(i);
            else if(!(g1->connect(curr_pair_index.first,nbor_pair_index.first)) || 
                    !(g2->connect(curr_pair_index.second,nbor_pair_index.second)))
                i = nbors.erase(i);
            else
                ++i;
        }
    } while(!nborsToBeGet.empty() && nbors.empty());

    if(nborsToBeGet.empty() && nbors.empty()) {
#ifdef DEBUG
        cout << "  ##no more neighbors##" << endl;
#endif
        if(!temp.empty()) {
#ifdef DEBUG
            cout << "  ##add path##" << endl;
            pair<vector<int>,vector<int> > p_indices = ag.toGraphIndices(temp);
            print_match(p_indices);
#endif
            result.push_back(temp);
        }
        return;
    }
    
    // split nbors into groups of possible mapping
    pair<vector<int>,vector<int> > nbor_indices = ag.toGraphIndices(nbors);
    vector<vector<int> > groups;
    bitset<NUMNBORS> notAssign;
    for(vector<int>::size_type i=0; i<nbors.size(); ++i)
        notAssign.set(i);
    while(notAssign.any()) {
        // get nbors with minimum graph index(1st)
        int min_index = 1000;
        for(vector<int>::size_type i=0; i<nbors.size(); ++i) {
            if(!notAssign.test(i))
                continue;
            if(nbor_indices.first[i] < min_index)
                min_index = nbor_indices.first[i];
        }
        vector<int> tmp_candidate;
        for(vector<int>::size_type i=0; i<nbors.size(); ++i) {
            if(!notAssign.test(i))
                continue;
            if(nbor_indices.first[i] == min_index) {
                tmp_candidate.push_back(nbors[i]);
                notAssign[i] = 0;
            }
        }
        extendGroup(ag,groups,tmp_candidate);
    }
    for(vector<vector<int> >::iterator i=groups.begin(); i!=groups.end(); ++i) {
#ifdef DEBUG
        cout << "  got one possible match extended from " << curr_id+1 << "(" << curr_pair_index.first+1 << "," << curr_pair_index.second+1 << ")" << endl;
        cout << "  AG:";
        for(vector<int>::iterator ii=i->begin(); ii!=i->end(); ++ii)
            cout << " " << *ii + 1;
        cout << endl;
        pair<vector<int>,vector<int> > pair_indices = ag.toGraphIndices(*i);
        print_match(pair_indices);
#endif
        for(vector<int>::iterator k=i->begin(); k!=i->end(); ++k) {
            pair<int,int> tmp_pair_index = ag.toGraphIndex(*k);
            visited1[tmp_pair_index.first] = true;
            visited2[tmp_pair_index.second] = true;
            nborsToBeGet.push_back(*k);
            temp.push_back(*k);
        }
        extend(ag,g1,g2,visited1,visited2,result,temp,nborsToBeGet);
        for(vector<int>::iterator k=i->begin(); k!=i->end(); ++k) {
            pair<int,int> tmp_pair_index = ag.toGraphIndex(*k);
            visited1[tmp_pair_index.first] = false;
            visited2[tmp_pair_index.second] = false;
            nborsToBeGet.pop_back();
            temp.pop_back();
        }
    }
    for(vector<int>::reverse_iterator ii=deleted_ids.rbegin(); ii!=deleted_ids.rend(); ++ii)
        nborsToBeGet.push_front(*ii);
}
// static method called by 'extendSCCS'
// make new path by concatenate each pair of pathes from 'final_result' and 'result'
static void addPathes(vector<vector<int> > &final_result, const vector<vector<int> > &result)
{
    if(result.empty())
        return ;
    for(vector<vector<int> >::const_iterator i=result.begin(); i!=result.end(); ++i)
        final_result.push_back(*i);
    /*
    vector<vector<int> > new_result;
    for(vector<vector<int> >::iterator i=final_result.begin(); i!=final_result.end(); ++i) {
        for(vector<vector<int> >::const_iterator j=result.begin(); j!=result.end(); ++j) {
            vector<int> temp(i->begin(),i->end());
            for(vector<int>::const_iterator k=j->begin(); k!=j->end(); ++k) {
                if(find(i->begin(),i->end(),*k) == i->end())
                    temp.push_back(*k);
            }
            new_result.push_back(temp);
        }
    }
    final_result.clear();
    for(vector<vector<int> >::iterator i=new_result.begin(); i!=new_result.end(); ++i)
        final_result.push_back(*i);
    */
}
static vector<vector<int> > extendSCCS(const AssocGraph &ag, 
        const MCSGraph *g1, const MCSGraph *g2,
        const pair<vector<int>,vector<int> > &sccs)
{
    int nodes1 = g1->numNodes();
    int nodes2 = g2->numNodes();
    vector<bool> visited1(nodes1,false);
    vector<bool> visited2(nodes2,false);
    vector<int> compsub;
    for(vector<int>::size_type i=0; i<sccs.first.size(); ++i) {
        compsub.push_back(ag.toAGIndex(sccs.first[i],sccs.second[i]));
        visited1[sccs.first[i]] = true;
        visited2[sccs.second[i]] = true;
    }
#ifdef DEBUG
    cout << endl << "in 'extendSCCS'" << endl;
    cout << "  AG: ";
    for(vector<int>::size_type i=0; i<compsub.size(); ++i)
        cout << " " << compsub[i]+1;
    cout << endl;
    print_match(sccs);
#endif
    
    vector<vector<int> > final_result;
    final_result.push_back(compsub);
    for(vector<int>::iterator i=compsub.begin(); i!=compsub.end(); ++i) {
#ifdef DEBUG
        pair<int,int> p_index = ag.toGraphIndex(*i);
        cout << "extending path from " << *i+1 << " (" << p_index.first+1 << ", " << p_index.second+1 << ")" << endl;
#endif
        deque<int> nborsToBeGet;
        nborsToBeGet.push_back(*i);
        vector<int> temp(compsub.begin(),compsub.end());
        vector<vector<int> > result;
        extend(ag,g1,g2,visited1,visited2,result,temp,nborsToBeGet);
        addPathes(final_result,result);
#ifdef DEBUG
        cout << "after extending and append to compsub" << endl;
        int num = 0;
        for(vector<vector<int> >::iterator ii=final_result.begin(); ii!=final_result.end(); ++ii) {
            pair<vector<int>,vector<int> > p_indices = ag.toGraphIndices(*ii);
            cout << "- path " << ++num << endl;
            print_match(p_indices);
        }
#endif
    }
    return final_result;
}

static bool inQuasiMCS(const vector<pair<vector<int>,vector<int> > > &quasiMCS,
        const vector<int> &g1, const vector<int> &g2)
{
    for(vector<pair<vector<int>,vector<int> > >::const_iterator i=quasiMCS.begin();
            i!=quasiMCS.end(); ++i) {
        if(i->first.size()==g1.size()) {
            bool same = true;
            for(vector<int>::size_type j = 0; j<g1.size(); ++j) {
                if(g1[j]!=i->first[j] || g2[j]!=i->second[j]) {
                    same = false;
                    break;
                }
            }
            if(same)
                return true;
        }
    }
    return false;
}

// remove duplicates from 'cliques'. called by 'findMCS'
static void removeDuplicate(vector<vector<int> > &cliques)
{
    // sort
    for(vector<vector<int> >::iterator i=cliques.begin(); i!=cliques.end(); ++i)
        sort(i->begin(),i->end());
    // remove duplicates
    vector<bool> remove(cliques.size(),false);
    for(vector<vector<int> >::size_type i=0; i<cliques.size(); ++i) {
        if(remove[i])
            continue;
        for(vector<vector<int> >::size_type j=i+1; j<cliques.size(); ++j) {
            if(remove[j])
                continue;
            if(cliques[i].size() != cliques[j].size())
                continue;
            bool same = true;
            for(vector<int>::size_type m=0; m<cliques[i].size(); ++m) {
                if(cliques[i][m] != cliques[j][m]) {
                    same = false;
                    break;
                }
            }
            if(same)
                remove[j] = true;
        }
    }
    int j = 0;
    for(vector<vector<int> >::iterator i=cliques.begin(); i!=cliques.end(); ) {
        if(remove[j])
            i = cliques.erase(i);
        else
            ++i;
        ++j;
    }
}

vector<pair<vector<int>,vector<int> > > findMCS(const MCSGraph *g1, const MCSGraph *g2,
        int rmax, int smin,
        float (*atomMapWeight)(const MCSGraph *g1, int i1, const MCSGraph *g2, int i2))
{
    // 1. AG and call clique finding algorithm
    AssocGraph ag(g1,g2,atomMapWeight);
    vector<vector<int> > cliques = findCliques(ag.getMatrix(),rmax);
    //cerr << "totally " << cliques.size() << " cliques" << endl;
    // 2. index of cliques -> index of g1&g2
    //    and remove SCCSs whose cardinality is lower than smin
    //    and retain unique SCCSs
    vector<pair<vector<int>,vector<int> > > quasiMCS;
    for(vector<vector<int> >::iterator i=cliques.begin(); i!=cliques.end(); ++i) {
        // graph indices of clique[i]
        pair<vector<int>,vector<int> > pair_index = ag.toGraphIndices(*i);
        // separate clique[i] into SCCSs
        vector<vector<int> > connect_index = connectGraph(g1,pair_index.first);
        // only retain those SCCSs with number of nodes >= smin
        for(vector<vector<int> >::iterator j=connect_index.begin(); j!=connect_index.end(); ++j) {
            if(static_cast<int>(j->size()) < smin)
                continue;
            vector<int> tmp_g1;
            vector<int> tmp_g2;
            for(vector<int>::iterator k=j->begin(); k!=j->end(); ++k) {
                tmp_g1.push_back(pair_index.first[*k]);
                tmp_g2.push_back(pair_index.second[*k]);
            }
            if(!inQuasiMCS(quasiMCS,tmp_g1,tmp_g2))
                quasiMCS.push_back(make_pair(tmp_g1,tmp_g2));
        }
    }
    //cerr << "after generating quasiMCS, there are " << quasiMCS.size() << " quasiMCSs" << endl;
    // 3. extend SCCSs
    cliques.clear();  // vector<vector<int> >
#ifdef DEBUG
    int ii = 0;
#endif
    for(vector<pair<vector<int>,vector<int> > >::iterator i=quasiMCS.begin(); i!=quasiMCS.end(); ++i) {
#ifdef DEBUG
        cout << "extend SCCS " << ++ii;
#endif
        vector<vector<int> > temp = extendSCCS(ag,g1,g2,*i);
        for(vector<vector<int> >::iterator j=temp.begin(); j!=temp.end(); ++j)
            cliques.push_back(*j);
#ifdef DEBUG
        cout << "==> totally " << temp.size() << " SCCSs" << endl << endl;
#endif
    }
    //cerr << "after extending SCCSs, there are " << cliques.size() << " cliques" << endl;
    // remove duplicates among 'cliques'
    removeDuplicate(cliques);
    // convert AG index to graph pair index
    quasiMCS.clear(); // vector<pair<vector<int>,vector<int> > >
    for(vector<vector<int> >::iterator i=cliques.begin(); i!=cliques.end(); ++i)
        quasiMCS.push_back(ag.toGraphIndices(*i));
    //cerr << "after converting index" << endl;
    // 4. calculate the weight of each mapping
    vector<float> sums;
    for(vector<pair<vector<int>,vector<int> > >::iterator i=quasiMCS.begin(); i!=quasiMCS.end(); ++i) {
        float sum = 0.;
        for(vector<int>::size_type j=0; j<(i->first).size(); ++j)
            sum += atomMapWeight(g1,i->first[j],g2,i->second[j]);
        sums.push_back(sum);
    }
    //cerr << "after calculating weight" << endl;
    // 5. extract those mapping(s) with maximum weight
    float max_sum = *std::max_element(sums.begin(),sums.end());
    vector<pair<vector<int>,vector<int> > > result;
    for(vector<pair<vector<int>,vector<int> > >::size_type i=0; i<quasiMCS.size(); ++i) {
        if(sums[i] == max_sum)
            result.push_back(quasiMCS[i]);
    }
    return result;
}

vector<int> mapWithMinChange(const MCSGraph *g1, const MCSGraph *g2,
        const vector<pair<vector<int>,vector<int> > > &mcs)
{
    // 1. minimize "number of mapping atoms that have non-mapping neighbors"
    vector<int> num_change(mcs.size());
    for(vector<pair<vector<int>,vector<int> > >::size_type i=0; i<mcs.size(); ++i) {
        int sum = 0;
        // to check neighbors of each mapping atom pair
        for(vector<int>::size_type j=0; j<mcs[i].first.size(); ++j) {
            vector<int> nbors_first = g1->getNborsOf(mcs[i].first[j]);
            bool miss = false;
            for(vector<int>::iterator k=nbors_first.begin(); k!=nbors_first.end(); ++k) {
                if(find(mcs[i].first.begin(),mcs[i].first.end(),*k) == mcs[i].first.end()) {
                    miss = true;
                    break;
                }
            }
            if(miss)
                ++sum;
            else {
                vector<int> nbors_second = g2->getNborsOf(mcs[i].second[j]);
                for(vector<int>::iterator k=nbors_second.begin(); k!=nbors_second.end(); ++k) {
                    if(find(mcs[i].second.begin(),mcs[i].second.end(),*k) == mcs[i].second.end()) {
                        miss = true;
                        break;
                    }
                }
                if(miss)
                    ++sum;
            }
        }
        num_change[i] = sum;
    }
    int min_change = *std::min_element(num_change.begin(),num_change.end());
    vector<int> min_index;
    for(vector<int>::size_type i=0; i<num_change.size(); ++i) {
        if(min_change == num_change[i])
            min_index.push_back(i);
    }
    // 2. minimize "number of bond changes within mapping atoms"
    vector<int> bond_change(min_index.size());
    for(vector<int>::size_type i=0; i<min_index.size(); ++i) {
        int sum = 0;
        for(vector<int>::size_type j=0; j<mcs[min_index[i]].first.size()-1; ++j) {
            for(vector<int>::size_type k=j+1; k<mcs[min_index[i]].first.size(); ++k) {
                int btype1 = (*g1)(mcs[min_index[i]].first[j],mcs[min_index[i]].first[k]);
                int btype2 = (*g2)(mcs[min_index[i]].second[j],mcs[min_index[i]].second[k]);
                if(btype1 && (btype1 != btype2))
                    ++sum;
            }
        }
        bond_change[i] = sum;
    }
    int min_bond_change = *std::min_element(bond_change.begin(),bond_change.end());
    vector<int> min_index1;
    for(vector<int>::size_type i=0; i<bond_change.size(); ++i) {
        if(min_bond_change == bond_change[i])
            min_index1.push_back(min_index[i]);
    }
    return min_index1;
}


