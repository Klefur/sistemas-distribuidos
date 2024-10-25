build: main.cpp
	g++ -fopenmp -o main.out main.cpp -lcfitsio
run: build
	./main.out -i elipsesimple1.fits -a 53 -r 50 -b 150 -u 6 -d 6
clean:
	rm -f *.out