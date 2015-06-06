#!/bin/bash
#Usage: ./rho Z N Emax J rhoFile
cd ../..
./rho 12 10 25 0 mg022rho.tsv
./rho 40 53 25 0 zr093rho.tsv
./rho 79 121 25 0 au200rho.tsv
mv mg022rho.tsv zr093rho.tsv au200rho.tsv scripts/rho
cd scripts/rho

echo "set terminal postscript enhanced dashed lw 1 'Helvetica' 14">gnuplotRho.gp
echo "set logscale y">>gnuplotRho.gp
echo "set xlabel 'E/MeV'">>gnuplotRho.gp
echo "set ylabel 'Arbitrary units'">>gnuplotRho.gp


echo "set output 'rho.ps'">>gnuplotRho.gp
echo "plot './mg022rho.tsv' u 1:2 w l, \\" >> gnuplotRho.gp

echo "'./zr093rho.tsv' u 1:2 w l, \\" >>gnuplotRho.gp

echo "'./au200rho.tsv' u 1:2 w l" >> gnuplotRho.gp 

gnuplot gnuplotRho.gp 