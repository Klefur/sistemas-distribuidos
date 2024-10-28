build: main.cpp
	g++ -fopenmp -o main.out main.cpp -lcfitsio
run: build
	./main.out -i elipsesimple1.fits -a 50 -r 40 -b 128 -u 4 -d 4
clean:
	rm -f *.out