# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= g++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG1 = main.x
PROG2 = functiontest.x


all: ${PROG1} ${PROG2}

clean:
			rm -f *.o

${PROG1} :	jacobiMethod.h jacobiMethod.cpp main.cpp
			${CPPflags} $^ -o ${PROG1} ${LIB}

${PROG2} : functiontest.cpp jacobiMethod.cpp jacobiMethod.h
			${CPPflags} $^ -o ${PROG2} ${LIB}
