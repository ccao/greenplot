include ../make.sys

objects = constants.o  hamk.o  input.o  para.o  wanndata.o

mod.a: $(objects)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

constants.o: constants.f90

para.o: para.f90 constants.o

input.o: input.f90 constants.o para.o

hamk.o: hamk.f90 constants.o input.o wanndata.o

wanndata.o: wanndata.f90 constants.o  para.o

clean:
	rm -rf *.x *.o *.mod *.a
