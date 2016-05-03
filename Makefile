MF=     Makefile
 
CC=     g++
 
CFLAGS= -g -fopenmp -D_USE_OMP -msse4.2 -O3 -fomit-frame-pointer -funroll-loops  
 
LFLAGS= -std=c++11 -I ./libFLASM -L ./libFLASM -Wl,-rpath=$(PWD)/libFLASM -lflasm -lahocorasick -I ./ahocorasick -L ./ahocorasick -I ./libdatrie/include -L ./libdatrie/lib -ldatrie -Wl,-rpath=$(PWD)/libdatrie/lib -lz -DNDEBUG -I ./libsdsl/include/ -L ./libsdsl/lib/ -lsdsl -ldivsufsort -ldivsufsort64 -Wl,-rpath=$(PWD)/libsdsl/lib
 
EXE=    bear
 
SRC=    bear.cc input.cc macsm.cc filter.cc aca.cc upgma.cc sw.cc nw.cc utils.cc matrices.cc csc.cc
 
HD=     beardefs.h globals.h filter.h aca.h EDNAFULL.h EBLOSUM62.h Makefile
 
# 
# No need to edit below this line 
# 
 
.SUFFIXES: 
.SUFFIXES: .cc .cpp .o 
 
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
	rm -r libdatrie
	rm -r libdatrie-0.2.8 
	rm -r libsdsl
	rm -r sdsl-lite
	rm -r libFLASM
