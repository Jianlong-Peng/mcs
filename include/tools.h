/*=============================================================================
#     FileName: tools.h
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Created: 2013-04-08 19:52:13
#   LastChange: 2013-04-11 10:37:17
#      History:
=============================================================================*/
#ifndef  MCS_TOOLS_H
#define  MCS_TOOLS_H

#include "graph.h"
#include <openbabel/atom.h>
#include <openbabel/mol.h>
int atomType(OpenBabel::OBAtom &atom, const char *atomtype_file);
// assign weight for node of association graph according to atom matching.
float weightOfNodeOfAG(const MCSGraph *g1, int i1, const MCSGraph *g2, int i2);

MCSGraph *OBMol2MCSGraph(OpenBabel::OBMol &mol, const char *atom_type_file);

#endif   /* ----- #ifndef MCS_TOOLS_H  ----- */

