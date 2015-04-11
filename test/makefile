CXX = g++
CFLAGS = -I ../include -Wall -O2

all: test_makeAG test_clique

test_makeAG: test_makeAG.cpp ../objs/graph.o
	$(CXX) $(CFLAGS) test_makeAG.cpp ../objs/graph.o -o test_makeAG

test_clique: test_clique.cpp ../objs/graph.o ../objs/clique.o
	$(CXX) $(CFLAGS) test_clique.cpp ../objs/graph.o ../objs/clique.o -o test_clique

