/*=============================================================================
#     FileName: clique.h
#         Desc: to find all maximal cliques of a given undirected graph
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-04-03 15:03:34
#   LastChange: 2013-04-09 22:08:09
#      History:
=============================================================================*/

#ifndef  MCS_CLIQUE_H
#define  MCS_CLIQUE_H

// ref: Bron C., Kerbosch J. Finding all cliques of an undirected graph. 1971
#include <vector>
#include "graph.h"
using std::vector;

// OBJ: to find all cliques for a given graph.
// ref: Born C. and Kerbosch J. Finding all cliques of an undirected graph. 1971
// Parameter
//   m   :
//   rmax: stop calculation of clique finding after 'rmax' recursion steps
//         default: 0, do not stop.
vector<vector<int> > findCliques(const MCSMatrix &m, int rmax=0);

// OBJ: to get atom mapping of the maximum clique
//      from all cliques, to get the one with maximum weight
// Parameters
//   - cliques    : result of 'findCliques'
//   - node_weight: used to calculate weight for each clique
//                  does here need 'node_weight' ??????
// Return
//   - AG indices of all clique(s) with maximum weight
vector<vector<int> > maxClique(const vector<vector<int> > &cliques,
        const vector<float> &node_weight);

#endif   /* ----- #ifndef MCS_CLIQUE_H  ----- */

