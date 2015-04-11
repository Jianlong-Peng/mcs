CXX = g++
OB_ROOT = /home/xmluo/jlpeng/program/openbabel
#CFLAGS = -I ./include -Wall -O2 -s
CFLAGS = -I ./include -Wall -g
LDFLAGS = -L $(OB_ROOT)/lib
OBJ_TARGET = objs
OBJS = $(OBJ_TARGET)/graph.o $(OBJ_TARGET)/clique.o $(OBJ_TARGET)/mcs.o $(OBJ_TARGET)/tools.o

all: $(OBJS) mcs mcs.bak

$(OBJ_TARGET)/graph.o: src/graph.cpp
	$(CXX) $(CFLAGS) -c src/graph.cpp -o $(OBJ_TARGET)/graph.o

$(OBJ_TARGET)/clique.o: src/clique.cpp
	$(CXX) $(CFLAGS) -c src/clique.cpp -o $(OBJ_TARGET)/clique.o

$(OBJ_TARGET)/mcs.o: src/mcs.cpp
	$(CXX) $(CFLAGS) -c src/mcs.cpp -o $(OBJ_TARGET)/mcs.o

$(OBJ_TARGET)/tools.o: src/tools.cpp
	$(CXX) $(CFLAGS) -I $(OB_ROOT)/include/openbabel-2.0 -c src/tools.cpp -o $(OBJ_TARGET)/tools.o

mcs: main.cpp $(OBJS)
	$(CXX) $(CFLAGS) -I $(OB_ROOT)/include/openbabel-2.0 $(LDFLAGS) $(OBJS) main.cpp -o mcs -lopenbabel

mcs.bak: main.bak.cpp $(OBJS)
	$(CXX) $(CFLAGS) -I $(OB_ROOT)/include/openbabel-2.0 $(LDFLAGS) $(OBJS) main.bak.cpp -o mcs.bak -lopenbabel

