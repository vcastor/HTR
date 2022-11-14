########################################################################
#              Makefile to compile Hartree-Fock-Roothaan
IDIR =../include
BDIR =..
ODIR = obj
SDIR =.
gg = gfortran -I$(IDIR)                                        #compiler
LIBS = -llapack -lblas                                        #libreries
FLAG = -fcheck=all -Wall

SRCF90 = $(wildcard *.f90)

OBJ = $(patsubst %.f90,$(ODIR)/%.o,$(SRCF90))

$(ODIR)/%.o: %.f90
	$(gg) -c -o $@ $< $(FLAG)

$(BDIR)/xRHFR.exe: $(OBJ)
	$(gg) -o $@ $^ $(LIBS) $(FLAG)
	@echo " ------------------------------------------------------ "
	@echo "   HartreeFockRoothaan has been successfully compiled"
	@echo " ------------------------------------------------------ "


clean:
	rm -f $(ODIR)/*.o $(BDIR)/*.exe
