MF=     Makefile
 
CC=     g++
 
CFLAGS= -g -fopenmp -D_USE_OMP -D_USE_MPFR -msse4.2 -O3 -fomit-frame-pointer -funroll-loops  
 
LFLAGS= -std=c++11 -lahocorasick -I ./ahocorasick -L ./ahocorasick -I ./libdatrie/include -L ./libdatrie/lib -ldatrie -Wl,-rpath=$(PWD)/libdatrie/lib -lz 
 
EXE=    bear
 
SRC=    bear.cc input.cc macsm.cc filter.cc aca.cc maxshift.cc upgma.cc sw.cc nw.cc utils.cc matrices.cc
 
HD=     beardefs.h globals.h filter.h aca.h EDNAFULL.h EBLOSUM62.h Makefile
 
# 
# No need to edit below this line 
# 
 
.SUFFIXES: 
.SUFFIXES: .cc .o 
 
OBJ=    $(SRC:.cc=.o) 
 
.cc.o: 
	$(CC) $(CFLAGS)-c $(LFLAGS) $< 
 
all:    $(EXE) 
 
$(EXE): $(OBJ) 
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS) 
 
$(OBJ): $(MF) $(HD) 
 
clean: 
	rm -f $(OBJ) $(EXE) *~

clean-all: 
	rm -f $(OBJ) $(EXE) *~
	rm -r ahocorasick
	rm -r libdatrie-0.2.8 
