<img width="200" src="conflow_logo_small.png" alt="ConFlow"/>
  
# ConFlow: Super Granular Convective Flow Generator 
  
[Predictive Science Inc.](https://www.predsci.com)  

--------------------------------  

## OVERVIEW ##
  
<img width="200" align=right src="doc/SG_sphere_lon0_lat60.png" alt="ConFlow"/>ConFlow computes a sequence of velocity maps for supergranule flows on the surface of the Sun.  Such photospheric velocity fields are essential for developing and testing realistic flux transport models and the analysis techniques for observational data.  

Conflow generates these maps analytically by specifying the convection spectrum of poloidal and toroidal modes (e.g., [Hathaway (1988)](https://doi.org/10.1007/BF00147251), [Hathaway et al. (2010)](https://doi.org/10.1088/0004-637X/725/1/1082)), and advecting them with the Sun’s axisymmetric differential rotation and meridional flows.

Conflow creates the velocity maps on a spherical surface grid of phi (longitude) and theta (colatitude). It can output the maps on a [staggered grid](https://github.com/predsci/ConFlow/blob/main/doc/psi_hipft_grid.png) in [HDF5](https://www.hdfgroup.org/solutions/hdf5) format or on an unstaggered grid in binary format. The staggered grid is designed to be directly used in  [OFT](https://github.com/predsci/oft)'s surface flux transport code [HipFT](https://github.com/predsci/hipft).  When using the maps with HipFT, one should set HipFT's meridional and differential flow coefficients to those used in the ConFlow computation.  

--------------------------------  
   
## HOW TO BUILD CONFLOW ##
  
ConFlow has been tested to work using GCC's `gfortran` (>8) compiler.  
It is recommended to use the latest compiler version available.  

HipFT requires the [HDF5](https://www.hdfgroup.org/solutions/hdf5) library.  
The library must be compiled by the same compiler ConFlow is using.  
When using the default system compiler (often GCC), the system's HDF5 dev packages should work.  

1. Find the build script from the `build_examples` folder that is closest to your setup and copy it into the top-level directory.  
2. Modify the script to set the `HDF5` library paths/flags and compiler flags compatible with your system environment.  
3. Modify the script to set the compiler options to reflect your setup.  
4. If using GCC's auto parallelization, set your `OMP_NUM_THREADS` environment variable to the number of threads you want to run with.  
5. Run the build script (for example, `./my_build.sh`).  
6. It is recommended to add the `bin` folder to your system path.  

--------------------------------  

## HOW TO RUN CONFLOW ##
  
### Setting Input Options  
  
`ConFlow` uses a namelist in an input text file.  
The name for the input text file must be set to `conflow.dat`  
  
A full working input file with all the default parameter options is provided in the file:  
  
`doc/conflow.dat.documentation`  
   
A detailed description of each parameter is also given in that file, and (in addition to this README) is the current main documentation of the code.  
  
We have also provided example input file for a use case in the `examples/` folder.  

### Launching the Code ###
    
To run `ConFlow`, set the desired run parameters into a file called  `conflow.dat`, then copy or link the `conflow` executable into the same directory as the input file and run the command:  
  
`./conflow`  
  
The code can be multi-threaded.  This requires proper setting of the `OMP_NUM_THREADS` environment variable before compilation when using GCC and its auto parallelization `ftree-parallelize-loops` flag.  

### Solution Output ###
  
The output of ConFlow (using a staggered grid) are HDF5 `vt` and `vp` velocity component map files in longitude-colatitude coordinates.  
A CSV file called `flow_output_list.csv` is generated that lists the output files with the time of output.  

--------------------------------




