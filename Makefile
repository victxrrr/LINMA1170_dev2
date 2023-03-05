CC=gcc
CFLAGS=-O2 -Wall

solve_deformation: matrix.c lu.c elasticity.c solve_deformation.c
	$(CC) $(CFLAGS) -o $@ $^ -lm ../gmsh-sdk/lib/libgmsh.so -Wl,-rpath,../gmsh-sdk/lib
	./solve_deformation square.geo 0.04
	rm -f solve_deformation

plot: plot.py matrix.csv solution.csv
	python3 $^

clean:
	rm -f solve_deformation
	rm -f *.o
	rm -f *.exe