SYSTEM = x86-64_linux
LIBFORMAT = static_pic

# ISIMA

# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------
CCC = g++
CC  = gcc
CCOPT = -fPIC -fexceptions -DNDEBUG -DIL_STD
COPT  = -fPIC
# ---------------------------------------------------------------------
# Cplex and Concert paths
# ---------------------------------------------------------------------
CONCERTVERSION = concert
CPLEXVERSION = CPLEX_Studio1210

CONCERTDIR = /opt/ibm/ILOG/$(CPLEXVERSION)/$(CONCERTVERSION)
CONCERTINCDIR = $(CONCERTDIR)/include/
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CPLEXDIR = /opt/ibm/ILOG/$(CPLEXVERSION)/cplex
CPLEXINCDIR = $(CPLEXDIR)/include/
CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

# ---------------------------------------------------------------------
# Flags
# ---------------------------------------------------------------------
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -lm -lpthread
CFLAGS  = $(COPT)  -I$(CPLEXINCDIR)
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CPLEXINCDIR)ilcplex -I$(CONCERTINCDIR)


LEMONINCDIR = /opt/lemon/include/
LEMONLIBDIR = /opt/lemon/lib/
LEMONCFLAGS = -I$(LEMONINCDIR)
LEMONCLNFLAGS = -L$(LEMONLIBDIR) -lemon


BOOSTINCDIR = /mnt/c/soft/boost_1_71_0/
BOOSTCFLAGS = -I$(BOOSTINCDIR)
# ---------------------------------------------------------------------
# Comands
# ---------------------------------------------------------------------
PRINTLN = echo

#---------------------------------------------------------
# .cpp Files
#---------------------------------------------------------
CPPFILES = main.cpp RSA.cpp cplexForm.cpp ShortestPath.cpp ExtendedGraph.cpp Slice.cpp Demand.cpp PhysicalLink.cpp Instance.cpp CSVReader.cpp input.cpp
#---------------------------------------------------------
# Files
#---------------------------------------------------------
all: main

main:
	$(CCC) -c -Wall -g -std=c++11 $(CCFLAGS) $(LEMONCFLAGS) $(BOOSTCFLAGS) $(CPPFILES)
	$(CCC) $(CCFLAGS) *.o -g -o exec $(CCLNFLAGS) $(LEMONCLNFLAGS)  -ldl
	rm -rf *.o *~ ^

clean:
	rm -rf *.o main

