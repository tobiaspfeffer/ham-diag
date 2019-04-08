CPP = g++ 
CFLAGS = -std=c++0x -O3
TARGET = HAM_diag
OBJ = main.o tree.o diagMC.o

#how to compile?
%.o: %.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

#how to link?
all: $(OBJ) 
	$(CPP) -o $(TARGET) $(OBJ) 

#dependencies
./main.o: ./main.cpp ./partition.hpp ./leaf.hpp ./tree.hpp ./diagMC.hpp ./diagMC.cpp
./tree.o: ./tree.cpp ./tree.hpp ./leaf.hpp ./partition.hpp
./diagMC.o: ./diagMC.cpp ./diagMC.hpp ./tree.hpp ./tree.cpp ./leaf.hpp

#routines
clean:
	\rm $(OBJ) *~ $(TARGET)
