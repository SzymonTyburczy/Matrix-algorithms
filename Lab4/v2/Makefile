CXX = g++
CXXFLAGS = -O3 -std=c++17 -march=native
EIGEN_INC = -I/usr/include/eigen3

TARGET = hmatrix
SRC = zad2.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(EIGEN_INC) $(SRC) -o $(TARGET)
	@echo "---------------------------------------------------"
	@echo "Kompilacja zakończona pomyślnie!"
	@echo "Aby uruchomić program, wpisz: ./$(TARGET)"
	@echo "---------------------------------------------------"

clean:
	rm -f $(TARGET)

.PHONY: all clean