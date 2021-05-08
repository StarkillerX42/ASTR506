
FC=gfortran

FLAGS=-O3 -g

OBJS=n_bodies.o integrators.o

all: init integrators_tests n_body_tests solsys

#setup:
#	python setup.py build_ext

init:
	$(FC) $(FLAGS) -c integrators.f90 n_bodies.f90 integrators_tests.f90 n_body_tests.f90 solsys.f90

py_ext:
	f2py -c -m fortpy integrators.f90 n_bodies.f90 --f90flags="$(FLAGS)" --quiet

integrators_tests:
	$(FC) $(FLAGS) -o integrators_tests integrators_tests.o $(OBJS)

n_body_tests:
	$(FC) $(FLAGS) -o n_body_tests n_body_tests.o $(OBJS)

solsys:
	$(FC) $(FLAGS) -o solsys solsys.o $(OBJS)

clean:
	rm solsys; rm integrators_tests; rm n_body_tests; rm *.mod; rm *.*o*

run:
	./integrators_tests
	time ./n_body_tests
	./solsys