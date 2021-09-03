
CPP := g++
CPPFLAGS := -Wall -g -std=c++11

DEPS := Reads.h Hash.h Graph.h

OBJ := Reads.o Hash.o Graph.o barcode-groups.o 

%.o: %.cpp $(DEPS)
	$(CPP) $(CPPFLAGS) -c $< -o $@

barcode-groups: $(OBJ)
	$(CPP) $(CPPFLAGS) $(OBJ) -o barcode-groups 

