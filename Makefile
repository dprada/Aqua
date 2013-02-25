## Checking compilers and libraries

# f2py command
F2PY=
# fortran compiler command
FCOMP=
# fortran compiler for f2py (not manually given)
FTYPE=
# lapack libraries
LAPACK_LIBS=
# options fortran compiler (-checkall, -fast, ..)
FOPTS=
# additional fortran flags
FFLAGS=

ifneq ($(comp),)
FCOMP=$(comp)
endif

ifneq ($(opt),)
FOPTS=$(opt)
endif


# Detecting platform (Linux,Mac,Cygwin)
PLATF_TYPE=$(shell uname -s)
# Detecting 32 or 64 bits
MACHINE_TYPE=$(shell uname -m)

ifeq ($(PLATF_TYPE),Darwin)
PLATF_TYPE=Mac
endif
ifeq ($(PLATF_TYPE),Linux)
PLATF_TYPE=Linux
endif

# Detecting f2py if "F2PY= "
ifeq ($(F2PY),)
F2PY_IN=$(shell which f2py 2>/dev/null)
F2PY2_IN=$(shell which f2py2 2>/dev/null)
ifneq ($(F2PY_IN),)
F2PY=f2py
endif
ifneq ($(F2PY2_IN),)
F2PY=f2py2
endif
endif

# Detecting compiler if "FCOMP= "
ifeq ($(FCOMP),)
IFORT_IN=$(shell which ifort 2>/dev/null)
GFORTRAN_IN=$(shell which gfortran 2>/dev/null)
PGF90_IN=$(shell which pgf95 2>/dev/null)
ifneq ($(PGF90_IN),)
FCOMP=pgf95
endif
ifneq ($(GFORTRAN_IN),)
FCOMP=gfortran
endif
ifneq ($(IFORT_IN),)
FCOMP=ifort
endif
endif

ifeq ($(FCOMP),gfortran)
FTYPE=gnu95
endif
ifeq ($(FCOMP),ifort)
ifeq ($(MACHINE_TYPE),x86_64)
FTYPE=intelem
else
FTYPE=intel
endif
endif
ifeq ($(FCOMP),pgf95)
FTYPE=pg
endif

# Configuring lapack libraries if "LAPACK_LIBS= "
PATH_LIBS=/usr/lib /usr/lib32 /usr/lib64 /lib /lib32 /lib64 
PATH_LOCAL_LIBS=$(subst :, ,$(LIBRARY_PATH))
PATH_LOCAL_LD_LIBS=$(subst :, ,$(LD_LIBRARY_PATH))

LAPACK_IN=$(shell find $(PATH_LIBS) $(PATH_LOCAL_LIBS) $(PATH_LOCAL_LD_LIBS) -name liblapack* 2>/dev/null )
MKL_IN=$(shell find $(PATH_LIBS) $(PATH_LOCAL_LIBS) $(PATH_LOCAL_LD_LIBS) -name libmkl* 2>/dev/null )

ifneq ($(LAPACK_IN),)
LAPACK_IN=1
endif
ifneq ($(MKL_IN),)
MKL_IN=1
endif

ifeq ($(FCOMP),ifort)
ifeq ($(MKL_IN),1)
ifeq ($(MACHINE_TYPE),x86_64)
LAPACK_LIBS= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_def -lpthread
else
LAPACK_LIBS= -lmkl_intel -lmkl_sequential -lmkl_core -lmkl_def -lpthread # ??
#LAPACK_LIBS= -lmkl_intel -lmkl_sequential -lmkl_core -lpthread -lm
endif
else
ifeq ($(LAPACK_IN),1)
LAPACK_LIBS= -llapack
endif
endif
endif

ifeq ($(FCOMP),gfortran)
LAPACK_LIBS= -llapack
endif

# Correction of libxdrfile name:
lxdrlib = libxdrfile.so.4.0.0
ifeq ($(PLATF_TYPE),Mac)
lxdrlib = libxdrfile.4.0.0.dylib
endif

# Checking everything right
ifeq ($(FCOMP),)
$(warning No Fortran compiler detected!)
endif
ifeq ($(F2PY),)
$(warning No F2PY command detected!)
endif
ifeq ($(LAPACK_LIBS),)
$(warning No LAPACK_LIBS detected!)
endif

# Pipe of output messages

PWD=$(shell pwd)
SOUT=1>> $(PWD)/INSTALL.log 2>> $(PWD)/INSTALL.log
CHECK=0

default: options prexdr io_formats/libxdrfile.so pref90 f90_libraries fin

options:
	@ echo "-----------------------------------------------"
	@ echo "# COMPILING AQUA in" $(PLATF_TYPE)"-"$(MACHINE_TYPE) "WITH:"
	@ echo "-----------------------------------------------"
	@ echo "Fortran Compiler:" $(FCOMP) "("$(FTYPE)")"
	@ echo "Fortran 2 Python:" $(F2PY)
	@ echo "Lapack libraries:" $(LAPACK_LIBS)
	@ echo "Fortran options :" $(FOPTS)
	@ if [ -e INSTALL.log ]; then rm INSTALL.log; fi

pref90:
	@ echo "--------------------------------------------"
	@ echo "# COMPILING FORTRAN LIBRARIES:"

prexdr:
	@ echo "--------------------------------------------"
	@ echo "# COMPILING XDR LIBRARY:"

fin:
	@ echo "--------------------------------------------"
	@ echo "# Check INSTALL.log for further info"
	@ echo "--------------------------------------------"

io_formats/libxdrfile.so: xdrfile-1.1.1.tar.gz
	@ echo '>>>>>> Installing the xdr library...' > INSTALL.log
	@ tar -zxvf xdrfile-1.1.1.tar.gz 1>/dev/null 2>/dev/null
	@ cd xdrfile-1.1.1/ ; ./configure --prefix=$(PWD)/xdrfiles --enable-fortran F77=$(FCOMP) --enable-shared $(SOUT)
	@ cd xdrfile-1.1.1/ ; make install $(SOUT)
	@ cd xdrfile-1.1.1/ ; make test $(SOUT)
	@ rm -r xdrfile-1.1.1
	@ if grep ': FAILED' INSTALL.log 1>/dev/null ; then echo '> Error: check the file INSTALL.log'; fi
	@ if ! grep ': FAILED' INSTALL.log 1>/dev/null ; then echo '> io_formats/libxdrfile.so ...   OK';\
	cp xdrfiles/lib/$(lxdrlib) io_formats/libxdrfile.so; rm -r xdrfiles; fi

f90_libraries: fort_general.so fort_enm.so fort_net.so fort_math.so fort_anal_trajs.so \
	fort_kin_anal.so io_formats/libdcdfile.so io_formats/libbinfile.so io_formats/libcell2box.so

fort_general.so: fort_general.f90
	@ echo '>>>>>> Compiling' $@ >> INSTALL.log
	@ $(F2PY) --opt=$(FOPTS) --f90flags=$(FFLAGS) --fcompiler=$(FTYPE) -c -m fort_general fort_general.f90 $(LAPACK_LIBS) $(SOUT)
	@ if [ ! -e $@ ]; then echo '> Error compiling' $@ ': check the file INSTALL.log'; fi
	@ if [ -e $@ ]; then echo '>' $@ '...   OK'; fi

water.so: water.f90
	@ echo '>>>>>> Compiling' $@ >> INSTALL.log
	@ $(F2PY) --opt=$(FOPTS) --f90flags=$(FFLAGS) --fcompiler=$(FTYPE) -c -m water water.f90 $(LAPACK_LIBS) $(SOUT)
	@ if [ ! -e $@ ]; then echo '> Error compiling' $@ ': check the file INSTALL.log'; fi
	@ if [ -e $@ ]; then echo '>' $@ '...   OK'; fi

fort_enm.so: fort_enm.f90
	@ echo '>>>>>> Compiling' $@ >> INSTALL.log
	@ $(F2PY) --opt=$(FOPTS) --f90flags=$(FFLAGS) --fcompiler=$(FTYPE) -c -m fort_enm fort_enm.f90 $(LAPACK_LIBS) $(SOUT)
	@ if [ ! -e $@ ]; then echo '> Error compiling' $@ ': check the file INSTALL.log'; fi
	@ if [ -e $@ ]; then echo '>' $@ '...   OK'; fi

fort_net.so: fort_net.f90
	@ echo '>>>>>> Compiling' $@ >> INSTALL.log
	@ $(F2PY) --opt=$(FOPTS) --f90flags=$(FFLAGS) --fcompiler=$(FTYPE) -c -m fort_net fort_net.f90 $(LAPACK_LIBS) $(SOUT)
	@ if [ ! -e $@ ]; then echo '> Error compiling' $@ ': check the file INSTALL.log'; fi
	@ if [ -e $@ ]; then echo '>' $@ '...   OK'; fi

fort_math.so: fort_math.f90
	@ echo '>>>>>> Compiling' $@ >> INSTALL.log
	@ $(F2PY) --opt=$(FOPTS) --f90flags=$(FFLAGS) --fcompiler=$(FTYPE) -c -m fort_math fort_math.f90 $(LAPACK_LIBS) $(SOUT)
	@ if [ ! -e $@ ]; then echo '> Error compiling' $@ ': check the file INSTALL.log'; fi
	@ if [ -e $@ ]; then echo '>' $@ '...   OK'; fi

fort_anal_trajs.so: fort_anal_trajs.f90
	@ echo '>>>>>> Compiling' $@ >> INSTALL.log
	@ $(F2PY) --opt=$(FOPTS) --f90flags=$(FFLAGS) --fcompiler=$(FTYPE) -c -m fort_anal_trajs fort_anal_trajs.f90 $(LAPACK_LIBS) $(SOUT)
	@ if [ ! -e $@ ]; then echo '> Error compiling' $@ ': check the file INSTALL.log'; fi
	@ if [ -e $@ ]; then echo '>' $@ '...   OK'; fi

fort_kin_anal.so: fort_kin_anal.f90
	@ echo '>>>>>> Compiling' $@ >> INSTALL.log
	@ $(F2PY) --opt=$(FOPTS) --f90flags=$(FFLAGS) --fcompiler=$(FTYPE) -c -m fort_kin_anal fort_kin_anal.f90 $(LAPACK_LIBS) $(SOUT)
	@ if [ ! -e $@ ]; then echo '> Error compiling' $@ ': check the file INSTALL.log'; fi
	@ if [ -e $@ ]; then echo '>' $@ '...   OK'; fi

fort_mss.so: fort_mss.f90
	@ echo '>>>>>> Compiling' $@ >> INSTALL.log
	@ $(F2PY) --opt=$(FOPTS) --f90flags=$(FFLAGS) --fcompiler=$(FTYPE) -c -m fort_mss fort_mss.f90 $(LAPACK_LIBS) $(SOUT)
	@ if [ ! -e $@ ]; then echo '> Error compiling' $@ ': check the file INSTALL.log'; fi
	@ if [ -e $@ ]; then echo '>' $@ '...   OK'; fi

io_formats/libdcdfile.so: io_formats/libdcdfile.f90
	@ echo '>>>>>> Compiling' $@ >> INSTALL.log
	@ cd $(PWD)/io_formats; $(F2PY) --opt=$(FOPTS) --f90flags=$(FFLAGS) --fcompiler=$(FTYPE) -c -m libdcdfile libdcdfile.f90 $(SOUT)
	@ if [ ! -e $@ ]; then echo '> Error compiling' $@ ': check the file INSTALL.log'; fi
	@ if [ -e $@ ]; then echo '>' $@ '...   OK'; fi

io_formats/libbinfile.so: io_formats/libbinfile.f90
	@ echo '>>>>>> Compiling' $@ >> INSTALL.log
	@ cd $(PWD)/io_formats; $(F2PY) --opt=$(FOPTS) --f90flags=$(FFLAGS) --fcompiler=$(FTYPE) -c -m libbinfile libbinfile.f90 $(SOUT)
	@ if [ ! -e $@ ]; then echo '> Error compiling' $@ ': check the file INSTALL.log'; fi
	@ if [ -e $@ ]; then echo '>' $@ '...   OK'; fi

io_formats/libcell2box.so: io_formats/libcell2box.f90
	@ echo '>>>>>> Compiling' $@ >> INSTALL.log
	@ cd $(PWD)/io_formats; $(F2PY) --opt=$(FOPTS) --f90flags=$(FFLAGS) --fcompiler=$(FTYPE) -c -m libcell2box libcell2box.f90 $(SOUT)
	@ if [ ! -e $@ ]; then echo '> Error compiling' $@ ': check the file INSTALL.log'; fi
	@ if [ -e $@ ]; then echo '>' $@ '...   OK'; fi

clean:
	@ rm *.so io_formats/*.so
	@ echo "Clean"
