all: euler solsys

#setup:
#	python setup.py build_ext

#test: class_01_28

euler:
	gfortran -O3 -o euler fort_integrators.f90

solsys:
	gfortran -O3 -o solsys solsys.f90
