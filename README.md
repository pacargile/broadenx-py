# broadenx-py
Python hooks and glue code for Bob Kurucz broadenx.for program

Files:

broaden.f90 -> Fortran code that is directly taken from Bob's broadenx.for code<br>
broaden.py  -> Python code using ctypes to access the hooks placed in the fortan<br>

Steps to install:

1) Compile the fortran code into a shared object library with something like:

gfortran -shared -fPIC -o libbroaden.so broaden.f90

2) Edit the path to this shared library within the broaden.py python code

Example of how to run the code within python:
  
  wl = np.array(wl) # a numpy array (or list) of wavelengths<br>
  flx = np.array(flx) # a numpy array (or list) of fluxes<br>
  res = 3000000.0 # initial resolution of spectrum<br>
  broaddict = {type:'GAUSSIAN','units':'KM','val':3.0} # info for a gaussian broadening of 3 km/s<br>
  
  # initialize the class <br>
  brd = broaden.broaden()<br>
  # do broadening and output a spectrum in a dictionary with keys "WAVE" & "FLUX"<br> 
  outspec = brd.broaden(wl,flx,res,broaddict)<br>
