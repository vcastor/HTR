########################################################################
#              Makefile to compile Hartree-Fock-Roothaan
gg = gfortran                                                  #compiler
lapack = -llapack                                         #library name
objects = HartreeFockRoothaan

$(objects): %: %.f90
	$(gg) -o $@.exe $< -llapack
	@echo " ------------------------------------------------------ "
	@echo "   HartreeFockRoothaan has been successfully compiled"
	@echo " "
	@echo "             run: ./launcher.sh input"
	@echo " ------------------------------------------------------ "


clean:
	rm -f $(objects).exe *.o *.mod *.a
