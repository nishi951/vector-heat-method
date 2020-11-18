# Vector Heat Method
## CS 468 Fall 2020 Final Project
N. Sharp, Y. Soliman, and K. Crane, “The vector heat method,” ACM Trans. Graph., vol. 38, no. 3, 2019. 

See the original paper [here.](https://www.cs.cmu.edu/~kmcrane/Projects/VectorHeatMethod/index.html)

## About
This repo contains a bare-bones julia implementation of the vector heat method for triangle meshes, created for the final project of CS 468 in Fall 2020.

You can run the code by modifying `vhm.jl` to specify a triangle mesh file (see
`meshes/` for some sample ones). A Julia environment is provided in the `vhm/`
folder; once initialized, running the main script is as easy as running

``` sh
julia --project=vhm -i vhm.jl
```

The `-i` flag will drop you into the Julia interpreter after the script
completes, keeping the plots open and allowing further investigation.

## API
The `preprocess()` function in `vhm_utils.jl` takes a list of vertices and triangle
faces, specified as `[i, j, k]` triples of vertex indices. After
checking for consistent orientation and removing duplicate vertices,
it computes all the relevant angles and tangent spaces needed to form the
laplacian and connection laplacian for triangle meshes.

After that, the system is formulated and solved in `vhm.jl` using Julia's
backslash operator.

Visualization helper functions are provided in `visualize_vhm.jl`.

`make_flatmesh.jl` creates a simple planar mesh useful for testing.



