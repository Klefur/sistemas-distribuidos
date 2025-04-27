build: main.cpp
	mpic++ -fopenmp -o main.out main.cpp 
run: build
	mpirun -np 4 ./main.out 100000 2
clean:
	rm -f *.out