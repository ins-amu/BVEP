BUILDDIRBLAS=.

SRCBLAS+=$(wildcard $(BUILDDIRBLAS)/*.f)

OBJFILESBLAS := $(SRCBLAS:.f=.o)
DEPFILESBLAS := $(SRCBLAS:.f=.d)

%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@ -lm

libblas.a : $(OBJFILESBLAS) 
	$(AR) rv libblas.a $^


clean :
	rm -f $(OBJFILESBLAS)
	rm -f libblas.a
