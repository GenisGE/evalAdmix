CXX=g++
CC=gcc

FLAGS=-O3

all : Cinterface

filereader_and_conversions.o : filereader_and_conversions.cpp filereader_and_conversions.h types.h
	$(CXX) -c filereader_and_conversions.cpp ${FLAGS}

extractors.o : extractors.cpp extractors.h types.h
	$(CXX) -c extractors.cpp ${FLAGS}

asort.o : asort.cpp asort.h types.h alloc.o
	$(CC) -c asort.cpp ${FLAGS}

alloc.o : alloc.h alloc.cpp types.h
	$(CC) -c -fPIC   alloc.cpp ${FLAGS}

evalAdmix.o : evalAdmix.h evalAdmix.cpp
	$(CC) -c evalAdmix.cpp ${FLAGS}

ngsevalAdmix.o : ngsevalAdmix.h ngsevalAdmix.cpp
	$(CC) -c ngsevalAdmix.cpp ${FLAGS}

Cinterface : Cinterface.cpp Cinterface.h filereader_and_conversions.o extractors.o asort.o alloc.o evalAdmix.o ngsevalAdmix.o
	$(CXX)  -o evalAdmix Cinterface.cpp  filereader_and_conversions.o extractors.o asort.o alloc.o evalAdmix.o ngsevalAdmix.o ${FLAGS} -lz -lpthread 

clean :
	rm -f *.o evalAdmix
