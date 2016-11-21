LIBNAME=repeated

OBJS= chidden.o cphidden.o cutil.o dist.o gar.o hidden.o kcountb.o \
	kserieb.o logitord.o binnest.o gausscop.o eigen.o \
	gaps.o calcs.o calcs1.o calcs2.o calcs3.o calcs4.o

$(LIBNAME)$(SHLIB_EXT):	$(OBJS)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) -o $@ $(OBJS) $(FLIBS)

clean:
	@rm -f *.o *.$(SHLIB_EXT)

realclean: clean
