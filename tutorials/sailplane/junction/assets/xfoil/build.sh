#!/bin/bash
# Build XFOIL 6.99 from Drela's source WITHOUT the floating-point traps of the
# Debian package (its gfortran build carries -ffpe-trap=invalid,zero and dies
# inside PACC with graphics off, which is how AeroSandbox's wrapper drives it).
# Double precision (-fdefault-real-8), X11 plot library as shipped (the
# wrapper switches graphics off with PLOP G; a display is never needed).
# Result: ./Xfoil/bin/xfoil -- point asb.XFoil(xfoil_command=...) at it.
set -euo pipefail
cd "$(dirname "$0")"
[ -d Xfoil ] || { [ -f xfoil6.99.tgz ] || curl -sSL -o xfoil6.99.tgz https://web.mit.edu/drela/Public/web/xfoil/xfoil6.99.tgz; tar xzf xfoil6.99.tgz; }
cd Xfoil/plotlib
cp config.make.gfortranDP config.make
make -s clean >/dev/null 2>&1 || true
make -s FFLAGS="-m64 -O2 -fdefault-real-8 -std=legacy -fallow-argument-mismatch" CFLAGS="-m64 -O2 -DUNDERSCORE -I/usr/include"
cd ../bin
make -s clean >/dev/null 2>&1 || true
make -s xfoil FC=gfortran CHK= DBL=-fdefault-real-8 \
     FFLAGS="-O2 -fdefault-real-8 -std=legacy -fallow-argument-mismatch" \
     FFLOPT="-O2 -fdefault-real-8 -std=legacy -fallow-argument-mismatch" \
     CFLAGS="-O2 -DUNDERSCORE" PLTLIB="-lX11" INSTALLCMD=true
ls -l "$PWD/xfoil"
