inverse: inverse.f90 cel.f r.f xigel.f90
	gfortran -g -Wall xigel.f90 inverse.f90 cel.f r.f -o inverse
