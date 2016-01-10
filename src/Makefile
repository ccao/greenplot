include ../make.sys

objects = calc_bulkgf.o calc_surfgf.o calc_gf.o

greenplot.x : greenplot.f90 $(objects)
	$(F90) -o greenplot.x greenplot.f90 $(INCLUDES) $(objects) ../modules/mod.a $(LAPACK95LIBS) $(BLAS95LIBS) $(MKLLIBS)

clean:
	rm -rf *.x *.o *.mod *.a
