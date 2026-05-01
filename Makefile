CXX = g++
CXXFLAGS = -O3 -fopenmp -std=c++17 -Iinclude

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

TARGET = mc

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f src/*.o $(TARGET)
