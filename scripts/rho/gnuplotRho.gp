set terminal postscript enhanced mono dashed lw 1 'Helvetica' 14
set logscale y
set xlabel 'E/MeV'
set ylabel 'Arbitrary units'
set output 'np238rho.ps'
plot './np238rho.tsv' u 1:2 w l
