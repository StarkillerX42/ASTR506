
FC=gfortran

FLAGS=-O3 -g

OBJS=integrators.o

all: init euler solsys

#setup:
#	python setup.py build_ext

#test: class_01_28

test:
	./integrators_tests; ./solsys

init:
	$(FC) $(FLAGS) -c integrators.f90 integrators_tests.f90 solsys.f90

euler:
	$(FC) $(FLAGS) -o integrators_tests integrators_tests.o $(OBJS)

solsys:
	$(FC) $(FLAGS) -o solsys solsys.o $(OBJS)

clean:
	rm euler; rm solsys; rm integrators_tests; rm *.mod; rm *.o; rm *.out
