## `OceanMesh2D: Precise distance-based two-dimensional automated mesh generation toolbox intended for coast ocean/shallow water flow models`

<p align="center">
  <img src = "nesting.png"> &nbsp &nbsp &nbsp &nbsp
  <img src = "ResoNa.png"> &nbsp &nbsp &nbsp &nbsp
  <img src = "Globalocean.jpg"> &nbsp &nbsp &nbsp &nbsp
</p>
OceanMesh2D is a set of user-friendly MATLAB functions to generate two-dimensional (2D) unstructured meshes for coastal ocean circulation problems. These meshes are based on a variety of feature driven geometric and bathymetric mesh size functions, which are generated according to user-defined parameters. Mesh generation is achieved through a force-balance algorithm combined with a number of topological improvement strategies aimed at improving the worst case triangle quality. The software embeds the mesh generation process into an object-orientated framework that contains pre- and post-processing workflows, which makes mesh generation flexible, reproducible, and script-able. 

## `Code framework` 
`OceanMesh2D`  consists of four standalone classes that are called in sequence. It requires no paid toolboxes to build meshes and has been tested to work with a trial version of MATLAB.

    OceanMesh2D::
    ├── geodata -- process geospatial data.
    ├── edgefx  -- build mesh size functions.
    ├── meshgen -- generate mesh based on mesh size functions and boundaries.
    └── msh     -- store, write, read, inspect, and visualize meshes and their axuillary components for numerical simulation.

## `Starting Out`

Clone or download and unzip the current <a href="https://github.com/CHLNDDEV/OceanMesh2D/archive/master.zip">repository</a>, 

Read the user guide available here:
https://www.overleaf.com/read/hsqjhvtbkgvj#/54715995/ (dynamic version, click download PDF), or https://doi.org/10.13140/RG.2.2.21840.61446/2 (static version)

The data for the following examples can be downloaded here: 
 https://drive.google.com/open?id=1LeQJFKaVCM2K59pKO9jDcB02yjTmJPmL
 
and here: http://www.soest.hawaii.edu/pwessel/gshhg/ (download the GSHH global coastline ESRI shapefile zip archive)
 
and here: ftp://topex.ucsd.edu/pub/srtm15_plus/ (download the SRTM15_PLUS global topobathy DEM "topo15_compressed.nc")
```
Featured in  ┌╼ Example_1_NZ.m   %<- A simple mesh around South Island New Zealand that uses GSHHS shoreline. 
user guide   ├── Example_2_NY.m   %<- A high-resolution mesh around the New York/Manhattan area that uses a DEM created from LiDAR data.  
             └── Example_3_ECGC.m %<- Builds a mesh for the western North Atlantic with a local high-resolution nest around New York
Featured in         ┌╼ Example_4_PRVI.m %<- Builds a mesh for the western North Atlantic with three high-resolution nests around Peurto Rico and US Virgin Islands
Geoscientific Model ├── Example_5_JBAY.m %<- An extremely high-fidelity (15-m) mesh from LiDAR data around Jamaica Bay with CFL-limiting.
Development paper[1]└── Example_6_GBAY.m %<- An example of the polyline/thalweg mesh size function along the Houston Ship Channel. 

```

## `References!`

If you make use of `OceanMesh2D` please include a reference to the following:
```

[1] - Roberts, K. J., Pringle, W. J., and Westerink, J. J., 2018. 
      OceanMesh2D 1.0: MATLAB-based software for two-dimensional unstructured mesh generation in coastal ocean modeling, 
      Geosci. Model Dev. Discuss., https://doi.org/10.5194/gmd-2018-203, in review.
[2] - Roberts, K. J., Pringle, W. J, 2018. 
      OceanMesh2D: User guide - Precise distance-based two-dimensional automated mesh generation toolbox intended for coastal
      ocean/shallow water. https://doi.org/10.13140/RG.2.2.21840.61446/2.       

```
