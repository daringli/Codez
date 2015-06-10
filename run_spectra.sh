datadir=/n/data2/bstefan
egdir=/n/home/bstefan/m/deexcite/final

plotdir=/n/home/bstefan/m/deexcite/final/scripts/plot_chart
#Z=10-90
#N=10-130
Z=10-90
N=10-130

for E in {5..20..5}
do
single_spectra=$E-talys-spectra.tsv

./nuclist -Z $Z -N $N -E $E -J 0 | ./spectra -t -l 10 -J 20 -d single > $single_spectra
cp $single_spectra $plotdir
done