SYSTEM = x86-64_linux
LIBFORMAT = static_pic

# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------
CCC = g++ -O0 -std=c++11
CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -Wno-ignored-attributes
# ---------------------------------------------------------------------
# Cplex, Concert, Lemon and Boost paths
# ---------------------------------------------------------------------
CONCERTVERSION = concert
CPLEXVERSION = CPLEX_Studio1210

CONCERTDIR = /home/mariabazotte/opt/ibm/ILOG/$(CPLEXVERSION)/$(CONCERTVERSION)
CONCERTINCDIR = $(CONCERTDIR)/include/
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CPLEXDIR = /home/mariabazotte/opt/ibm/ILOG/$(CPLEXVERSION)/cplex
CPLEXINCDIR = $(CPLEXDIR)/include/
CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

LEMONINCDIR = /home/mariabazotte/lemon-1.3.1/include/
LEMONLIBDIR = /home/mariabazotte/lemon-1.3.1/lib/
LEMONCFLAGS = -I$(LEMONINCDIR)
LEMONCLNFLAGS = -L$(LEMONLIBDIR) -lemon

BOOSTINCDIR = /home/mariabazotte/boost_1_71_0/
BOOSTCFLAGS = -I$(BOOSTINCDIR)


# ---------------------------------------------------------------------
# Flags
# ---------------------------------------------------------------------
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread -ldl
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -lm -lpthread
CFLAGS  = $(COPT)  -I$(CPLEXINCDIR)
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)

#---------------------------------------------------------
# .cpp Files
#---------------------------------------------------------
#CPPFILES = main.cpp formulation/*.cpp solver/*.cpp tools/*.cpp topology/*.cpp
CPPFILES = main2.cpp lagrangian/*.cpp lagrangian/solver/*.cpp lagrangian/formulation/*.cpp  formulation/*.cpp solver/*.cpp tools/*.cpp topology/*.cpp subgradient.cpp
# ---------------------------------------------------------------------
# Comands
# ---------------------------------------------------------------------
PRINTLN = echo

#---------------------------------------------------------
# Files
#---------------------------------------------------------
all: main

main:
	doxygen ../dconfig
	$(CCC) -c -Wall -g $(CCFLAGS) $(LEMONCFLAGS) $(BOOSTCFLAGS) $(CPPFILES)
	$(CCC) $(CCFLAGS) *.o -g -o exec $(CCLNFLAGS) $(LEMONCLNFLAGS)
	rm -rf *.o *~ ^

clean:
	rm -rf *.o main ../Output/LP/* ../Output/*.csv ../doc/* out

