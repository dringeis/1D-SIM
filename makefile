## INSTRUCTION: make EXE='executable_filename' ice

# Makefile for my sea-ice model
# compilation by 
#FC	=	f77 # fortran 77
FC	=	gfortran #ifc -w # fortran 95

# -w above means that the warning messages are ignored...there are a few left
# but they don't affect the compilation

#F	=	ifc -r8 -O3 -tpp7 -W0 -c  test done with Jan...not working
#LIBDIR = -L/home/hdx1/gavin/lib
LIBS  = #-llapack -lblas #-lnag
#FFLAGS 	= -s -O -Pv -lc    # used preprocessing
#FFLAGS  = -g  #  used for debug
#FFLAGS  = -pg  #  used for profiling run the code and then type gprof
#FFlAGS  = -O0  

EXE = 'zoupa'
SRC	= src/ice.f90 src/parameter_mod.f90 src/constants_mod.f90 src/global_var_mod.f90 src/option_mod.f90 src/ini_get.f90 src/wind_forcing.f90 src/ice_strength.f90 src/bvect.f90 src/viscouscoefficient.f90 src/Cw_coefficient.f90 src/Fu.f90 src/SOR.f90 src/upwind_adv.f90 src/util.f90 src/prep_fgmres_NK.f90 src/identity.f90 src/JacfreeVec.f90 src/fgmresD.f90 src/dcopy.f90 src/ddot.f90 src/daxpy.f90 src/EVP2solver.f90 src/EVP2Dprecond.f90 src/output_results.f90 src/output_file.f90 src/sigma_s.f90 src/Fu_EVP.f90


ice: $(SRC) $(COMMON) 
	$(FC) $(FFLAGS) -o $(EXE) $(SRC)  #$(LIBS)


#lib: ice
#	ar rv libice.a $(OBJ)

tidy:	
	rm ./*.o src/*~ ./zoupa


# $(LIBDIR)  $(LIBS)

 
#$@
