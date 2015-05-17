To obtain and use these code perform the following steps:

```
mkdir sbp_projection
cd sbp_projection
git clone https://github.com/bfam/sbp_projection_operators
git clone https://github.com/bfam/sbp_projection_2d
cd sbp_projection_2d
```
To use the DG with Hesthaven and Warburton's codes (known
to work with commit 42b078b3b0) complete the following
```
cd sbp_projection
git clone git@github.com:tcew/nodal-dg.git
cd sbp_projection_2d/src
cp ../../nodal-dg/Codes1.1/Codes2D/MakeCylinder2D.m MakeCurved2D.m
patch MakeCurved2D.m < MakeCurved2D.patch
tar xjf straight_v3.tar.bz
cd ..
```

To run the codes launch MATLAB from the `sbp_projection_2d` directory and run

```matlab
addpath drivers
```

The SBP driver scripts in the 'drivers' directory have the following naming
convention
```
driver_{CONFIG}_p{ORDER}_{INTERFACE}.m
```
where {CONFIG} is a '2block' or '3block' configuration, {ORDER} is the interior
SBP finite difference order, {INTERFACE} 'c', 'n', and 'u' result in
a conforming, nested, or unnested interface, respectively.

The SBP-DG driver scripts in the 'drivers' directory have the following naming
convention
```
driver_sbpdg_p{SBPORDER}_p{DGORDER}.m
```
where {SBPORDER} is the interior SBP finite difference order and {DGORDER} is
the DG polynomial order.

To generate latex tables and tikz figures with error and convergence rates run
```
make_tables
make_plots
```

To generate an eigenvalue spectrum for a coupled SBP-DG simulation with both an
upwinded and central flux run
```
driver_make_eigenvalues
```

References:
```
@BOOK{HesthavenWarburton2008,
  title = {Nodal Discontinuous {G}alerkin Methods: {A}lgorithms, Analysis, and
           Applications},
  publisher = {Springer},
  year = {2008},
  author = {Hesthaven, Jan S. and Warburton, Tim},
  volume = {54},
  series = {Texts in Applied Mathematics},
  doi = {10.1007/978-0-387-72067-8}
}
```
