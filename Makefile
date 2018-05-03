CXX=mpicc

CXXFLAGS=-std=c++11

TARGET=distributed_mat_mul

distributed_mat_mul:distributed_mat_mul.o

$(TARGET): $(TARGET).cpp
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(TARGET).o

clean:
	rm *.o distributed_mat_mul 
