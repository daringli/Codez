datadir=/n/data2/bstefan
egdir=/n/home/bstefan/m/deexcite/git/11-05-2015

plotdir=/n/home/bstefan/m/deexcite/git/11-05-2015/scripts/plot_bars
#Z=10-90
#N=10-130
Z=5-20
A=28

for E in {20..20..1}
do
part_spectra=$E-t-spectra.tsv

./nuclist -Z $Z -A $A -E $E -J 0 | ./spectra -l 10 -J 20 -t -d particles > $part_spectra
cp $part_spectra $plotdir
done