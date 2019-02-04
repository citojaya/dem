# Project: mill
# Makefile created by Dev-C++ 4.9.9.2

gcc   = g++
WINDRES = windres.exe
RES  = 
OBJ  =  Cell.o WornHoleOneCode.o Force.o  InputData.o Mill.o Point.o Simulate.o $(RES)
LINKOBJ  = Cell.o WornHoleOneCode.o Force.o InputData.o Mill.o Point.o Simulate.o$(RES)
LIBS =  -lm -L"lib"   
INCS =  -I"include" 
CXXINCS =  -I"include" 
BIN  = mill 
CXXFLAGS = $(CXXINCS) -g3  
CFLAGS = -O -systype bsd43 $(INCS)  -g3
RM = rm -f

.PHONY: all all-before all-after clean clean-custom
all-before:
all-after:
clean-custom:

all: all-before mill all-after


clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(gcc) $(LINKOBJ) -o "mill" 

WornHoleOneCode.o: WornHoleOneCode.cpp
	$(gcc) -c WornHoleOneCode.cpp -o WornHoleOneCode.o $(CXXFLAGS)

Force.o: Force.cpp
	$(gcc) -c Force.cpp -o Force.o $(CXXFLAGS)

InputData.o: InputData.cpp
	$(gcc) -c InputData.cpp -o InputData.o $(CXXFLAGS)

Mill.o: Mill.cpp 
	$(gcc) -c Mill.cpp -o Mill.o $(CXXFLAGS)

Point.o: Point.cpp
	$(gcc) -c Point.cpp -o Point.o $(CXXFLAGS)

Simulate.o: Simulate.cpp
	$(gcc) -c Simulate.cpp -o Simulate.o $(CXXFLAGS)
Cell.o: Cell.cpp
	$(gcc) -c Cell.cpp -o Cell.o $(CXXFLAGS)
