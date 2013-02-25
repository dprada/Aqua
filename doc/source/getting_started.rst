.. _getting_started_pyn:

***************
Getting Started
***************

The whole project is public available at https://github.com/dprada/Pynoramix.git .
There is not stable version yet, for this reason the use of these libraries it is under your responsability.

Getting Pynoramix
=================

The testing version
+++++++++++++++++++

This version of pynoramix is the seed of the first stable version, and
thereby the most robust configuration at the moment.
All functions and attributes found on this documentation are available here.

The source code can be downloaded in the web page https://github.com/dprada/Pynoramix.git as a zip or tar.gz file.
Both options can be found with the options "zip" or "Downloads":

.. image:: _static/screenshot_github.png


The last developing version
+++++++++++++++++++++++++++

In this version new functions or corrections are updated almost
weekly. Because of this, keeping an eye on the history of the project
is highly recommended:

https://github.com/dprada/Pynoramix/commits/master

The last **developing** version of the project can be cloned with git:

.. sourcecode:: bash

   git clone -b master git://github.com/dprada/Pynoramix.git

The use of git is recommended since libraries can be easily updated.

.. Todo:: In the future the project will be included in easy_install
   or setup.py
   (http://packages.python.org/an_example_pypi_project/setuptools.html#registering-your-project)




Installing
===========

Pynoramix depends on some packages:

- Python 2.7
- Fortran Compiler (gfortran, intel fortran compiler, ...)
- NumPy
- f2py
- python-dev (to fix the problem with file Python.h)
- liblapack and liblapack-dev (or similar: atlas, blas, mkl, ...)
- Pylab (recommended but not necessary)

Mac users:

- Install gfortran via this link (http://gcc.gnu.org/wiki/GFortranBinaries#MacOS)
- Install f2py via macports: sudo port install py-f2py
- Remember to use python2.7 instead of the stock python



After solving the dependencies, the Makefile needs to be executed to compile the fortran core of Pynoramix.
This installation script has some variables which can be optionally fullfilled manually:

.. sourcecode:: bash

   F2PY=             # f2py command (f2py,f2py2,...)
   FCOMP=            # fortran compiler command (gfortran, ifort,...)
   FTYPE=            # fortran compiler for f2py (not manually given)
   LAPACK_LIBS=      # lapack libraries (-llapack, -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_def -lpthread, ...)
   FOPTS=            # options of the fortran compiler used (-fast, -checkall, ...)
   FFLAGS=           # additional fortran flags

If these variables are left in blank, they will be detected automatically. 
At this point, and in the directory of Pynoramix, the following command needs to be executed:

.. sourcecode:: bash

   make

If the installation run without troubles, Pynoramix is ready to be used.

.. warning:: Do not forget to add Pynoramix to your python path: export PYTHONPATH=$PYTHONPATH:/path/to/Pynoramix


Being updated
=============

The last modifications of the developing version can be easily downloaded with git.
The command 'git pull' executed over the Pynoramix directory checks and obtains the changes.
Once this has been done, compiling the changed libraries is mandatory. 
Since the Makefile script detects the modified files, running it again is enough.

.. sourcecode:: bash

   git pull
   make


If at any moment the installation needs to be done from scratch, the
following command removes the compiled files:

.. sourcecode:: bash

   make clean

