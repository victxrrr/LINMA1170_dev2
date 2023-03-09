CC=gcc
CFLAGS=-O2 -Wall

solve_deformation: matrix.c lu.c elasticity.c solve_deformation.c
	$(CC) $(CFLAGS) -o $@ $^ -lm ../gmsh-sdk/lib/libgmsh.so -Wl,-rpath,../gmsh-sdk/lib
	./solve_deformation square.geo 0.04
	rm -f solve_deformation

plot: plot.py permuted_matrix.csv bandK.csv
	python3 $^

plot_perm: plot.py permuted_matrix.csv permuted_solution.csv
	python3 $^

plot_K_permK: plot.py matrix.csv permuted_matrix.csv
	python3 $^

plot_Ksol_permKsol: plot.py solution.csv permuted_solution.csv
	python3 $^

main: matrix.c lu.c
	$(CC) $(CFLAGS) -o $@ $^ -lm
	./main

clean:
	rm -f solve_deformation
	rm -f solve_deformation_band
	rm -f *.o
	rm -f *.exe