########################################################################
#              Makefile to compile Hartree-Fock-Roothaan
IDIR =../include
BDIR =../bin
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
	@echo " "
	@echo "             run: ./launcher.sh input"
	@echo " ------------------------------------------------------ "


clean:
	rm -f $(OBDIR).exe *.o *.mod *.a