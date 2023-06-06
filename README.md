# MPI EXERCISE 2

Finite difference methods 1D Poisson Implementation and 2D Poisson Implementation via MPI

## How to use
### 0. Load MPI module
$module load cports openmpi

### 1. Build the project
$ make clean && make

### 2. Run Q.1-3 script on Chuck(1D Implementation)
$ mpirun -np [Number of processor] ./main [Number of Square Grids]
$ example: mpirun -np 4 ./main 13

### 3. Run Q.4 script on Chuck(2D Implementation)
$ mpirun -np [Number of Processor] ./main2d [Number of Square Grids]
$ example: mpirun -np 4 ./main2d 13

### 4. Change the grid size in "Poisson1d.h"

### 5. Repeat step.1 to step.3








