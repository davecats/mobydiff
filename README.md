# mobydiff
<img align="left" src="https://github.com/davecats/mobydiff/blob/main/.mobydiff.png" width="400"> <br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/>
A WIP simple and (reasonably) efficient solver of the incompressible Navier-Stokes equations

## Dependencies

*  _cmake > 3.22_
*  _gcc >= 13_

For CPU support
* _fftw >= 3.1_

For GPU support
* NVIDIA: NVHPC + CUDA _12+_

## Compile

Compile both CPU/GPU paths via 
```
./compile.sh all
```
