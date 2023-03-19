CC=gcc
CFLAGS=-O2 -Wall

solve_deformation: matrix.c lu.c elasticity.c solve_deformation.c
	$(CC) $(CFLAGS) -o $@ $^ -lm -Wno-unused-function ../gmsh-sdk/lib/libgmsh.so -Wl,-rpath,../gmsh-sdk/lib
	./solve_deformation square.geo 0.08
	rm -f solve_deformation

measurements: matrix.c lu.c elasticity.c solve_deformation.c
	for number in 1.0 0.5 0.3 0.2 0.1 0.09 0.08 0.05 0.04 0.03 0.025; do \
		$(CC) $(CFLAGS) -o $@ $^ -lm -Wno-unused-function ../gmsh-sdk/lib/libgmsh.so -Wl,-rpath,../gmsh-sdk/lib ; \
		./measurements square.geo $$number ; \
		rm -f measurements ; \
	done

plot: plot.py files/K228.csv files/lu_K228.csv files/bandK228.csv files/lu_bandK228.csv files/symK228.csv files/cholesky_symK228.csv
	python3 $^

time_plot: time_plot.py files/time_results.csv
	python3 $^

clean:
	rm -f solve_deformation
	rm -f *.o
	rm -f *.exe