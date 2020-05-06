gfortran -fPIC -O3 -cpp -DWITHOUTMPI -c -o cooling_module.o cooling_module.f90
gfortran -fPIC -O3 -cpp -DWITHOUTMPI -c -o cooling_ramses.o cooling_ramses.f90
f2py -m cooling cooling_ramses.f90 -h cooling.pyf --overwrite-signature
f2py -lgfortran -c *.o cooling.pyf