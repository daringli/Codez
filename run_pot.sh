atadir=/n/data2/bstefan
egdir=/n/home/bstefan/m/deexcite/git/11-05-2015

plotdir=/n/home/bstefan/m/deexcite/git/11-05-2015/scripts/potential
#Z=10-90
#N=10-130
Z=5-20
A=28

E=1
potfile=pots.tsv
transfile=trans.tsv

./nuclist -Z $Z -A $A -E $E -J 0 | ./pot -l 0 -e p -e n -e alpha > $potfile
./nuclist -Z $Z -A $A -E $E -J 0 | ./trans -l 0 -e p -e n -e alpha > $transfile
./nuclist -Z $Z -A $A -E $E -J 0 | ./sep -p p -p n -p alpha 

cp $potfile $transfile $plotdir
