#include<iostream>
#include<fstream>

using namespace std;

//TO INSTALL IN DEBIAN DO THIS:
//TO SUCCESSFULLY COMPILE YOU NEED GCC AND OPENMP INSTALLED
// sudo apt-get install gcc libgomp1

int main()
{
    ofstream file("makefile.mak");
    file << "CPP = c++" << endl;
    file << "CC = gcc" << endl;
    file << "RES = " << endl;
    file << "VPATH = src: src/core/: src/slvrs: " << endl;
    file << "OBJDIR = obj" << endl;

    //--------------------------------------------------------------------------
    // CORE FILES
    //--------------------------------------------------------------------------
    file << "OBJ =  $(OBJDIR)/ps_basic.o ";
    file << "$(OBJDIR)/ps_gridgen.o ";
    file << "$(OBJDIR)/ps_read.o ";
    file << "$(OBJDIR)/ps_implicit.o ";
    file << "$(OBJDIR)/ps_prntmesh.o ";
    file << "$(OBJDIR)/ps_prntpnts.o ";
    file << "$(OBJDIR)/mathops.o ";
    file << "$(OBJDIR)/settings.o ";
    file << "$(OBJDIR)/bc.o ";
    file << "$(OBJDIR)/krnl_readwrite.o ";
    file << "$(OBJDIR)/krnl_sprt.o ";
    file << "$(OBJDIR)/krnl_cmn.o ";
    file << "$(OBJDIR)/krnl_wls.o ";
    file << "$(OBJDIR)/krnl_mls.o ";
    file << "$(OBJDIR)/krnl_rbf.o ";
    file << "$(OBJDIR)/krnl_sph.o ";
    file << "$(OBJDIR)/krnl_gfd.o ";
    file << "$(OBJDIR)/intrf_basic.o ";
    file << "$(OBJDIR)/intrf_vel.o ";
    file << "$(OBJDIR)/intrf.o ";
    file << "$(OBJDIR)/trnsprt.o ";
    file << "$(OBJDIR)/ns.o ";
	file << "$(OBJDIR)/prtcl_shft.o ";

    //--------------------------------------------------------------------------
    // SOLVERS
    //--------------------------------------------------------------------------
    file << "$(OBJDIR)/slv_trnsprt.o ";
    file << "$(OBJDIR)/slv_ls.o ";
    file << "$(OBJDIR)/slv_pf.o ";
    file << "$(OBJDIR)/slv_ns1p.o ";
    file << "$(OBJDIR)/slv_ns2p.o ";
    file << "$(OBJDIR)/slv_nslgr.o ";
	
    //--------------------------------------------------------------------------
    // MAIN FILE
    //--------------------------------------------------------------------------
    file << "$(OBJDIR)/main.o " << " $(RES)" << endl;

    file << "LINKOBJ = $(OBJ)" << endl;
    file << "BIN  = run" << endl;
    file << "CXXFLAGS = $(CXXINCS)" << endl;
    file <<"CFLAGS = $(INCS)" << endl;

    file <<"RM = rm -f" << endl;

    file <<".PHONY: all clean clean-custom" << endl;

    file <<"all: run" << endl;

    file <<"clean: clean-custom" << endl;
    file <<"	${RM} $(OBJ) $(BIN)" << endl;

    file <<"$(BIN): $(OBJ)" << endl;
    file <<"	$(CPP) $(LINKOBJ) -o \"run\" $(LIBS) -std=c++17 -fopenmp -Wall -O3" << endl;

    file << "$(OBJDIR)/%.o: %.cpp" << endl;
    file <<"	$(CPP) -c $(CXXFLAGS) $< -o $@ -std=c++17 -fopenmp -Wall -O3" << endl;
    
    file.close();
    return 0;
}
