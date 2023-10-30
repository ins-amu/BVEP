LIB=portP/chk
FFLAGS=-O
OBJ=\
err1.o\
err2.o\
err3.o\
errk.o\
extr.o\
gama.o\
gmad.o\
lgac.o\
lgad.o\
lgea.o\
llza.o\
lrad.o\
lrpa.o\
ltqd.o\
ltsqc.o\
lyad.o\
lyma.o\
lyrc.o\
lysc.o\
mach.o\
mfed.o\
mfte.o\
mntb.o\
plya.o\
plyc.o\
prac.o\
prad.o\
prsa.o\
qbgd.o\
qblg.o\
qggd.o\
qgsg.o\
qpad.o\
qpra.o\
rnrm.o\
sdba.o\
sple.o\
stk1.o\
stk2.o\
svad.o\
svbd.o\
svcd.o\
svda.o\
svdb.o\
svdc.o\
thyp.o\
vdad.o\
vdbd.o\
vded.o

lib.a:	$(OBJ)
	update lib.a $?

install:	lib.a
	mv lib.a /usr/local/lib/lib$(LIB).a
	rm *.o

test: test.o
	f77 test.o -l$(LIB)
	time a.out
