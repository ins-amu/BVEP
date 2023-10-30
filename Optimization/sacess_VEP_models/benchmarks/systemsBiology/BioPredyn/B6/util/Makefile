# gen_deviates
GDOBJ = error.o mathLib.o ioTools.o global.o#distributions.o random.o dSFMT.o dSFMT_str_state.o 

# header files

HEADS = global.h error.h 
#LSA_HEADS = global.h sa.h MPI.h error.h 
#RND_HEADS = global.h error.h #random.h 
#DIS_HEADS = global.h error.h #distributions.h random.h 

#targets

all:  $(GDOBJ)
#deviates.o distributions.o ioTools.o error.o lsa.o random.o 

#gen_deviates: $(GDOBJ)
#	$(CC) -o gen_deviates $(CFLAGS) $(LDFLAGS) $(GDOBJ) $(LIBS)

#deviates.o:  $(DIS_HEADS) deviates.c 
#	$(CC) $(CFLAGS) -c deviates.c -o deviates.o

#distributions.o: $(DIS_HEADS) distributions.c
#	$(CC) $(CFLAGS) -c distributions.c -o distributions.o

ioTools.o: $(HEADS) ioTools.h ioTools.c
	$(CC) -L$(GSL)/lib $(CFLAGS) -c ioTools.c -o ioTools.o

mathLib.o: $(HEADS) mathLib.h mathLib.c
	$(CC) -L$(GSL)/lib $(CFLAGS) -c mathLib.c -o mathLib.o

error.o: $(HEADS) error.c
	$(CC) -L$(GSL)/lib  $(INCLUDES) $(CFLAGS) -c error.c -o error.o

global.o: $(HEADS) global.c
	$(CC) -L$(GSL)/lib  $(INCLUDES) $(CFLAGS) -c global.c -o global.o

#lsa.o: $(LSA_HEADS) lsa.c
#	$(CC) $(CFLAGS) -c lsa.c -o lsa.o

#random.o: $(RND_HEADS) random.c
#	$(CC) $(CFLAGS) -c random.c -o random.o
	
#dSFMT.o: $(RND_HEADS) dSFMT.c
#	$(CC) $(CFLAGS) -c dSFMT.c -o dSFMT.o	
	
#dSFMT_str_state.o: $(RND_HEADS) dSFMT_str_state.c
#	$(CC) $(CFLAGS) -c dSFMT_str_state.c -o dSFMT_str_state.o		

# parallel stuff

# lsa-mpi.o: lsa.c
#	$(MPICC) -c -o lsa-mpi.o $(MPIFLAGS) $(CFLAGS) lsa.c

# ... and here are the cleanup and make deps rules

clean:
	rm -f *.o core*

