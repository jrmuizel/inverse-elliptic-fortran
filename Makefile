inverse: inverse.f90 cel.f r.f
	gfortran -Wall xigel.f90 inverse.f90 cel.f r.f -o inverse
