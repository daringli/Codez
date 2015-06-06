datadir=/n/data2/bstefan
egdir=/n/home/bstefan/m/deexcite/git/11-05-2015

rootdir=/n/data2/bstefan
#Z=10-90
#N=10-130
Z=12
A=28

E=30
spectra_file=spectra.root

./nuclist -Z $Z -A $A -E $E -J 0 | ./spectra -t -l 10 -J 20 -d full > spectra-full.tsv
./out-spectra2root -f $spectra_file < spectra-full.tsv
cp $spectra_file $rootdir
