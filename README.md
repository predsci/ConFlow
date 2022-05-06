######![ConFlow](conflow_logo_small.png)
  
# ConFlow: Super Granular Convective Flow Generator 
  
## OVERVIEW ##
  
ConFlow creates 1024 by 512 velocity maps of supergranules.  
  
It constructs a spectrum of spherical harmonic amplitudes  
s(l,m) and calculates the components of the velocity field at  
an array of points in theta=colatitude and phi=longitude.  
  
The velocity vector (u,v) is in the (phi,theta) direction.  
  
The data arrays are written to direct access disk files.   
  
Differential rotation is simulated by evolving  
the spectral coefficients using 3-coefficient fits (North-South symmetric).  
  
Meridional flow is simulated by evolving  
the spectral coefficients using 6-coefficient fits.  
