build: main.cpp
	g++ -fopenmp -o main.out main.cpp -lcfitsio
run: build
	./main.out AS209_edge.fits
clean:
	rm -f *.out