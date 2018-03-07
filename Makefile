# -*-Makefile-*-

OBJS = main.o runProgram.o preprocessing.o queryValidation.o fileIO.o subsampling2.o correlation.o csdFiltering.o utils.o readFileTestUtils.o

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
	CC =/usr/local/opt/llvm/bin/clang++ -O3 -fopenmp -v
	DEBUG = -g
	CFLAGS = -Wall -c $(DEBUG)
	LFLAGS = -Wall $(DUBUG) -L /usr/local/opt/llvm/lib
else
	CC = g++ -O3 -fopenmp
	DEBUG = -g
	CFLAGS = -Wall -c $(DEBUG)
	LFLAGS = -Wall $(DUBUG)
endif


DiCoN: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o CSD-CS

main.o: main.cpp runProgram.h queryValidation.h
	$(CC) $(CFLAGS) main.cpp

preprocessing.o: preprocessing.h preprocessing.cpp utils.h myStructs.h
	$(CC) $(CFLAGS) preprocessing.cpp

queryValidation.o: queryValidation.h queryValidation.cpp myStructs.h
	$(CC) $(CFLAGS) queryValidation.cpp

runProgram.o: runProgram.h runProgram.cpp queryValidation.h preprocessing.h fileIO.h subsampling2.h correlation.h csdFiltering.h runProgram.h myStructs.h utils.h
	$(CC) $(CFLAGS) runProgram.cpp

fileIO.o: fileIO.h fileIO.cpp preprocessing.h readFileTestUtils.h utils.h myStructs.h
	$(CC) $(CFLAGS) fileIO.cpp

subsampling2.o: subsampling2.h subsampling2.cpp
	$(CC) $(CFLAGS) subsampling2.cpp

correlation.o: correlation.h correlation.cpp
	$(CC) $(CFLAGS) correlation.cpp

csdFiltering.o: csdFiltering.h csdFiltering.cpp myStructs.h
	$(CC) $(CFLAGS) csdFiltering.cpp

utils.o: utils.h utils.cpp myStructs.h
	$(CC) $(CFLAGS) utils.cpp

readFileTestUtils.o: readFileTestUtils.h readFileTestUtils.cpp
	$(CC) $(CFLAGS) readFileTestUtils.cpp

clean:
	\rm *.o CSD-CS
