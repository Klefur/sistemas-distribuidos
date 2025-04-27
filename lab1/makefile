build: main.cpp
	g++ -fopenmp -o main.out main.cpp -lcfitsio
run:
	./main.out -i elipsesimple1.fits -a 60 -r .4 -b 128 -u 1 -d 8
clean:
	rm -f *.out