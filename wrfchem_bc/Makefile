.SUFFIXES: .o .F .f .f90

F90    = ifort -g 
LIBS   = -L/opt/apps/intel13/netcdf/3.6.3/x86_64/lib -lnetcdf 
INCLUDE_MODULES = -I/opt/apps/intel13/netcdf/3.6.3/x86_64/include
FFLAGS =

EXEC = wrfchembc 


OBJS = module_GEOS_lib.o \
       module_wrfchem_lib.o \
       main_bc_wrfchem.o

##dependencies
#$(OBJECTS) : makefile


${EXEC} :       ${OBJS}
		${F90} -o $@ ${FFLAGS} ${OBJS} ${LIBS} 

.f.o:
		${F90} ${FFLAGS} -c $<
.f90.o:
		${F90} ${FFLAGS} -c ${INCLUDE_MODULES} $<

clean:
		rm -f core ${EXEC} ${OBJS} *.mod
