CC=gcc
CFLAGS=-O2 -Wall

solve_deformation: matrix.c lu.c elasticity.c solve_deformation.c
	$(CC) $(CFLAGS) -o $@ $^ -lm -Wno-unused-function ../gmsh-sdk/lib/libgmsh.so -Wl,-rpath,../gmsh-sdk/lib
	./solve_deformation square.geo 0.03
	rm -f solve_deformation

measurements: matrix.c lu.c elasticity.c solve_deformation.c
	for number in 1.0 0.9 0.7 0.5 0.3 0.1 0.05 0.04 0.03 0.025; do \
		$(CC) $(CFLAGS) -o $@ $^ -lm -Wno-unused-function ../gmsh-sdk/lib/libgmsh.so -Wl,-rpath,../gmsh-sdk/lib ; \
		./measurements square.geo $$number ; \
		rm -f measurements ; \
	done


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