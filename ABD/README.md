# Affine Body PBD

An example of how to run affine body pbd is in `pbd.m`.

Affine box box collison is done in `AffineBoxBox.cpp` and a simple test is in `testAffineBoxBox.m`.



## Dependencies
[Eigen 3.4.0](https://eigen.tuxfamily.org/index.php?title=Main_Page)  
[FCL](https://github.com/flexible-collision-library/fcl)  

## How to build AffineBoxBox.cpp
```
 mex COPTIMFLAGS='-O3 -DNDEBUG' -I"D:\ZiyanXiong\Cloth_Sim_Env\eigen-3.4.0" -I"C:\Program Files (x86)\libccd\include" -I"C:\Program Files (x86)\fcl\include" -L"C:\Program Files (x86)\libccd\lib" -L"C:\Program Files (x86)\fcl\lib" AffineBoxBox.cpp fcl.lib ccd.lib
```

Potential errors and solutions:

1. `MATLAB curshes after running testAffineBoxBox.m where there are contacts`

   There is a `ccd.dll` file in MATLAB's `bin` folder. But what we need is the `ccd.dll` in `libccd\bin` folder. Temporary solution: delete the `ccd.dll` file in MATLAB's `bin` folder and add the path to `libccd\bin` in system `PATH`.