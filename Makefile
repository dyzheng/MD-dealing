CPLUSPLUS=g++
OPTION=-O3

OBJS=main.o\
input.o\
matrix3.o\
gfun.o\
cell.o\
cellFile.o\
cellVASP.o\
cellLAMMPS.o\
cellABACUS.o\
atoms.o\
ext.o\
pdf.o\
pdf2d.o\
msd.o\
ssf.o\
dsf.o\
vel.o\
velcor.o\
powers.o\
vacuum.o\
data3D.o\
math.o\
pseudo.o\
iprof.o\
isf.o\
bdf.o\
average.o\
insert.o\

.cpp.o:
	${CPLUSPLUS} ${OPTION} -c $< -o $@


D310.exe :${OBJS}
	${CPLUSPLUS} ${OPTION} -o ABACUS_ANALYSE.exe ${OBJS}

clean:
	rm -f *.o *.exe
