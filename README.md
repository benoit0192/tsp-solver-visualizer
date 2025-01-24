# Travelling Salesman Problem (TSP) Solver & Visualizer
<div align="center">
  <img src="/asset/trefoil_knot.gif" alt="Trefoil Knot">
</div>

A simple TSP solver based on the Christofides algorithm is provided in this repository.<br>
It also includes a map visualizer implemented using OpenGL 3.3.<br>

To ensure reproducibility, a Docker environment setup is included.<br>
Using Docker is optional, as long as the necessary package dependencies are installed on your system.<br>

## Docker environment setup (optional)
### Build the docker image
```shell
$ docker/build.sh
```
### Enter the docker container
```shell
$ script/docker.sh
```
## Build the solver
Once built, the program `bin/vanilla-tsp' should be generated.<br>
```shell
$ make solver
```
## Run the solver
```shell
$ bin/vanilla-tsp <DATA_FILE> <OUT_DIR> --verbose
```
- `<DATA_FILE>` is the path to the data file containing the user-defined node data.
Please refer to the `examples` folder for the file format.
- `<OUT_DIR>` is optional and represents the folder directory where the solver will save the output data.
- `--verbose` is optional and will print extra informational messages.

## Build the map visualizer
Once built, the program `bin/map-renderer' should be generated.<br>
```shell
$ make map-renderer
```
## Run the map visualizer
```shell
$ bin/map-renderer <DATA_FILE>
```
- `<DATA_FILE>` is the path to the data file generated by the solver, containing the nodes and path information.
