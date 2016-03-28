# pruebo a hacer mi makefile
BIN_PATH=~/bin/tools/bin/
FC=ifort
OBJS1=rmsd.o pdb_in_out.o
OBJS2=rmsd.o pdb_in_out.o para_dihedros.o
LIBS=-L/opt/intel/mkl/10.0.1.014/lib/em64t -lmkl_lapack -lmkl -lguide -lpthread

default: rmsd_fit move_dihedral

rmsd_fit: main_rmsd.f90 $(OBJS1)
	@ $(FC) $(OBJS1) main_rmsd.f90 -o $(BIN_PATH)$@ $(LIBS)
	@echo "compilado rmsd_fit"

move_dihedral: move_dihedral.f90 $(OBJS2)
	@ $(FC) $(OBJS2) move_dihedral.f90 -o $(BIN_PATH)$@ $(LIBS)
	@echo "compilado move_dihedral.f90"

rmsd.o: rmsd.f90
	@ $(FC) -c rmsd.f90
	@echo "compilado el objeto rmsd.o"

pdb_in_out.o: pdb_in_out.f90
	@ $(FC) -c pdb_in_out.f90
	@echo "compilado el objeto pdb_in_out.o"

para_dihedros.o: para_dihedros.f90
	@ $(FC) -c para_dihedros.f90
	@echo "compilado el objeto para_dihedros.o"

clean: 
	@rm *.o *.mod $(BIN_PATH)rmsd_fit $(BIN_PATH)move_dihedral
	@echo "Limpio"

ayuda:
	@echo "Compila rmsd_fit junto con los objetos rmsd.o y pdb_in_out.o"

purgar:
	@rm *~
	@echo "Limpio"