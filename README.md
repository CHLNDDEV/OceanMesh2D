<p align="center">
  <a href="https://github.com/CHLNDDEV/OceanMesh2D"><img alt="OceanMesh2D" src="imgs/Oceanmesh_logo_white.png" width="95%"></a>
  <p align="center">Precise distance-based two-dimensional automated mesh generation toolbox intended for coast ocean/shallow water flow models.</p>
</p>

Table of contents
=================

<!--ts-->
   * [OceanMesh2D](#oceanmesh2d)
   * [Table of contents](#table-of-contents)
   * [Getting help](#getting-help)
   * [Contributing](#contributing)
   * [Code framework](#code-framework)
   * [Starting out](#starting-out)
   * [Testing](#testing)
   * [References](#references)
   * [Gallery](#gallery)
   * [Changelog](#changelog)
<!--te-->

## IMPORTANT NOTE:
This is the default and recommended `PROJECTION` branch. Please use it unless you otherwise require legacy (`MASTER` branch) or the absolute newest features (`DEV` branch).

OceanMesh2D is a set of user-friendly MATLAB functions to generate two-dimensional (2D) unstructured meshes for coastal ocean circulation problems. These meshes are based on a variety of feature driven geometric and bathymetric mesh size functions, which are generated according to user-defined parameters. Mesh generation is achieved through a force-balance algorithm combined with a number of topological improvement strategies aimed at improving the worst case triangle quality. The software embeds the mesh generation process into an object-orientated framework that contains pre- and post-processing workflows, which makes mesh generation flexible, reproducible, and script-able.

Getting help
==============

Besides posting [issues](https://github.com/CHLNDDEV/OceanMesh2D/issues) with the code on Github, you can also ask questions via our Slack channel [here](https://join.slack.com/t/oceanmesh2d/shared_invite/zt-ldexcc7m-LkkxsUTUekHbJNxwlyrMYQ).

Note: If the slack link invite isn't working, please send either one of us and an email and we'll fix it. By default, the invite link expires every 30 days. 

Otherwise please reach out to either Dr. William Pringle (wpringle@anl.gov) or Dr. Keith Roberts (krober@usp.br) with questions or concerns or feel free to start an Issue in the issues tab above.


Contributing
============

All contributions are welcome!

To contribute to the software:

1. [Fork](https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github/fork-a-repo) the repository.
2. Clone the forked repository, add your contributions and push the changes to your fork.
3. Create a [Pull request](https://github.com/CHLNDDEV/OceanMesh2D/pulls)

Before creating the pull request, make sure that the examples pass...

Some things that will increase the chance that your pull request is accepted:
- Write minimal working [examples](https://en.wikipedia.org/wiki/Minimal_working_example#:~:text=In%20computing%2C%20a%20minimal%20working,to%20be%20demonstrated%20and%20reproduced) that demonstrate the functionality.
- Write good commit and pull request messages.

Code framework
================
`OceanMesh2D`  consists of four standalone classes that are called in sequence. It requires no paid toolboxes to build meshes and has been tested to work with a trial version of MATLAB.

    OceanMesh2D::
    ├── geodata -- process geospatial data.
    ├── edgefx  -- build mesh size functions.
    ├── meshgen -- generate mesh based on mesh size functions and boundaries.
    └── msh     -- store, write, read, inspect, and visualize meshes and their axuillary components for numerical simulation.

Starting Out
============

Clone or download and unzip the current [repository](https://github.com/CHLNDDEV/OceanMesh2D/archive/Projection.zip)

PLEASE READ THE USER GUIDE!
A recent pdf of the user guide is located in this branch. For a continually updated version click [here](https://www.overleaf.com/read/hsqjhvtbkgvj#/54715995/) (wait for compilation and then click download PDF)

Run the "setup.sh" bash script to download the required m_map package and the base datasets:
- GSHHG global shoreline 
- SRTM15_PLUS global topobathy DEM

Additional data required for some of the following examples must be downloaded manually from [here](https://drive.google.com/open?id=1LeQJFKaVCM2K59pKO9jDcB02yjTmJPmL). Specifically, Examples 2, 3, 4, 5 and 5b require additional datasets from the google drive folder while base datasets are sufficient for the other examples. 
```
Featured in  ┌╼ Examples/Example_1_NZ.m   %<- A simple mesh around South Island New Zealand that uses GSHHS shoreline.
user guide   ├── Examples/Example_2_NY.m   %<- A high-resolution mesh around the New York/Manhattan area that uses a DEM created from LiDAR data.
             └── Examples/Example_3_ECGC.m %<- Builds a mesh for the western North Atlantic with a local high-resolution nest around New York
Featured in         ┌╼ Examples/Example_4_PRVI.m %<- Builds a mesh for the western North Atlantic with three high-resolution nests around Peurto Rico and US Virgin Islands
Geoscientific Model ├── Examples/Example_5_JBAY.m %<- An extremely high-fidelity (15-m) mesh from LiDAR data around Jamaica Bay with CFL-limiting.
Development paper[1]└── Examples/Example_6_GBAY.m %<- An example of the polyline/thalweg mesh size function along the Houston Ship Channel.

```

See #testing to test OceanMesh2D on your system.

Testing
==========

To ensure the software is fully functional on your system before building some crazy meshes it is strongly recommended to run the tests (RunTests.m) in the Tests/ directory.

We test all pull requests using this test suite on a local host before accepting. For subtantial pull requests we will also test the Examples from the Examples/ directory. 

References!
==============

If you make use of `OceanMesh2D` please include a reference to [1], and to any of [2]-[5] if pertinent ([latex .bib file](https://github.com/CHLNDDEV/OceanMesh2D/tree/Projection/UserGuide/OceanMesh2D_library.bib)). We would also appreciate using our [logo](https://github.com/CHLNDDEV/OceanMesh2D/tree/Projection/imgs) in a presentation featuring `OceanMesh2D`.
```

[1] - Roberts, K. J., Pringle, W. J., and Westerink, J. J., 2019.
      OceanMesh2D 1.0: MATLAB-based software for two-dimensional unstructured mesh generation in coastal ocean modeling,
      Geoscientific Model Development, 12, 1847-1868. https://doi.org/10.5194/gmd-12-1847-2019.
[2] - Roberts, K. J., Pringle, W. J, 2018.
      OceanMesh2D: User guide - Precise distance-based two-dimensional automated mesh generation toolbox intended for coastal
      ocean/shallow water. https://doi.org/10.13140/RG.2.2.21840.61446/2.
[3] - Roberts, Keith J. Unstructured Mesh Generation and Dynamic Load Balancing for Coastal Ocean Hydrodynamic Simulation, 2019.
      PhD Thesis, University of Notre Dame. https://curate.nd.edu/show/4q77fr0022c.
[4] - Roberts, Keith J., Pringle W.J., Westerink J. J. Contreras, M.T., Wirasaet, D., 2019.
      On the automatic and a priori design of unstructured mesh resolution for coastal ocean circulation models,
      Ocean Modelling, 144, 101509. https://doi.org/10.1016/j.ocemod.2019.101509.
[5] - Pringle, W. J., Wirasaet, D., Roberts, K. J., and Westerink, J. J., 2020.
      Global Storm Tide Modeling with ADCIRC v55: Unstructured Mesh Design and Performance,
      Geoscientific Model Development Discussions. https://doi.org/10.5194/gmd-2020-123.

```
In addition, best practice when using software in a scientific publication is to cite the permanent doi corresponding to the version used (e.g., for reproducibility). All our releases are archived at the following `Zenodo` repository doi [link](https://doi.org/10.5281/zenodo.1341384).
```
Authors (202X). CHLNDDEV/OceanMesh2D: OceanMesh2D VX.X. Zenodo. https://doi.org/10.5281/zenodo.1341384
```
Please fill in the version (VX.X), author list and year corresponding to the version used.

## `DISCLAIMER: `
The boundary of the meshing domain must be a polygon (first point equals the last and non-self intersecting) but it does not need to be simplified. Read the user guide for more information about the inputs.

GALLERY:
=========
* These images can be made by running the [examples](https://github.com/CHLNDDEV/OceanMesh2D/tree/Projection/Examples)
<p align="center">
  <img src = "imgs/VizGilgo.png"> &nbsp &nbsp &nbsp &nbsp
  <img src = "imgs/Example1.png"> &nbsp &nbsp &nbsp &nbsp
  <img src = "imgs/nesting.png"> &nbsp &nbsp &nbsp &nbsp
  <img src = "imgs/ResoNA.png"> &nbsp &nbsp &nbsp &nbsp
  <img src = "imgs/Globalocean.jpg"> &nbsp &nbsp &nbsp &nbsp
</p>

Changelog
=========

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)

### Unreleased

## Added
- `mesh2d` interface improvements to filter small polygons.
- Support for creation of `fort.20` files for forcing rivers by @Jiangchao3
- Mesh "cleaning" modes moderate and aggressive transfer nodal attributes via improvements to `msh.map_mesh_properties`

## Changed
- `msh.plot()` overhaul. All options specified via kwarg.
- `msh.plot()` option `subset` option is now called `subdomain`
- `msh.plot()` arbitary f13 option now utilizes `colormap` kwarg
- `utilities/extract_subdomain` now is called with kwargs.
- Cleaning up `msh.bound_courant_number()` to use `msh.get_boundary_of_mesh()` for Delaunay-triangulation and allowing `msh.clean()` to do the transfer of attributes automatically.

## Fixed
- Boundary labeling fix
- Prompt when labeling bcs using `outer` kwarg in `make_bc`

## Deleted
- Deprecating `msh.CheckTimestep()` for `msh.bound_courant_number`. Added error message and instruction in the CheckTimestep help

### [3.3.0] - 2020-12-21
## Fixed
- Users without mapping toolbox could not read in shapefiles because of a bug that
  made them required to have a 3d shapefiles.
- plotting `gdat` with no shoreline.
- plotting a mesh's bathymetry with a non-zero datum using cmocean.
- cell-averaging interpolation method in msh.interp fixed for unequal lon-lat DEM grid spacings

## Added
- Mesh patch smoother
- Ability to remesh abritary patches of elements within the domain while respecting user-defined mesh sizes and the patches boundaries.
- Ability to use the TPXO9 Atlas for the tidal bcs and sponge (inside tidal_data_to_ob.m and Calc_Sponge.m) by using '**' wildcards in place of the constituent name within the tidal atlas filename (the atlas has an individual file for each constituent).
- Introducing 'auto_outer' option for the make_bc msh method which populates the bc for the outermost mesh boundary polygon (ignores islands)
- Changelog to README
- "mapMeshProperties" msh method ports over mesh properties for a mesh subset
- 'invert' option in the msh.interp method to turn off the DEM value inversion typically performed

## Changed
- for the make_bc msh method 'auto'/'auto_outer' options, allowing for the 'depth' method of classification to use the interpolated depths on the mesh if gdat is empty.
- improving help for make_bc msh method, Make_f15.m and Calc_Sponge.m
- renamed "ExtractSubDomain.m" to "extract_subdomain.m"
- improving "extract_subdomain.m" help and facilitating NaN-delimited polygons
- ability to return boundary as a cell in "get_boundary_of_mesh" msh method
- "Example_1_NZ.m" includes example of plotting bcs of a msh subset
- using "mapMeshproperties" method in "fixmeshandcarry"
- using "fixmeshandcarry" in the "cat" msh method
- improving warning and error messages for the "interp" msh method
- adding geofactor into "writefort15" for the GAHM vortex model
