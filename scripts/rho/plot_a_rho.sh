#!/bin/bash
#Usage: ./plotRho Z N J Emax rhoFile
cd ../..
file=np238rho.tsv
outfile=np238rho.ps

./rho 93 145 6 0 $file
mv $file scripts/rho
cd scripts/rho

echo "set terminal postscript enhanced mono dashed lw 1 'Helvetica' 14">gnuplotRho.gp
echo "set logscale y">>gnuplotRho.gp
echo "set xlabel 'E/MeV'">>gnuplotRho.gp
echo "set ylabel 'Arbitrary units'">>gnuplotRho.gp


echo "set output '$outfile'">>gnuplotRho.gp
echo "plot './$file' u 1:2 w l" >> gnuplotRho.gp


gnuplot gnuplotRho.gp 