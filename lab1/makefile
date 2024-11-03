build: main.cpp
	g++ -fopenmp -o main.out main.cpp -lcfitsio
run: build
	./main.out -i elipsesimple1.fits -a 50 -r .4 -b 128 -u 4 -d 1
run2: build
	./main.out -i AS209_edge.fits -a 300 -r .4 -b 1500 -u 4 -d 1
clean:
	rm -f *.out