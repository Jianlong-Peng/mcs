#include <iostream>
#include <vector>
#include "graph.h"

using namespace std;

int main()
{
	int labels[] = {1,1,1,1,1,1};
	MCSGraph g1(labels,2);
	MCSGraph g2(labels,3);
	
	// add edge
	// g1: 0-1
	// g2: 0-1-2
	g1(0,1) = 1;
	g2(0,1) = 1;
	g2(1,2) = 1;
	// display g1 and g2
	cout << "g1" << endl;
	g1.display(cout," ");
	cout << "g2" << endl;
	g2.display(cout," ");
	
	// association graph
	AssocGraph ag(&g1,&g2);
	cout << "\nassociation graph" << endl;
	ag.getMatrix().display(cout," ");
	cout << endl;
	
	return 0;
}


