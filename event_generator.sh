#! /bin/bash
#Run the gun!
datadir=/n/data2/bstefan

egdir=/n/home/bstefan/m/deexcite/git/11-05-2015
quasi_file=qout.tsv
output=dqout.tsv
gun_file=test.gun
root_file=rootfile.root
fire_file=fired.root
events=200
gg_events=4200

GG=/net/home/bstefan/ggland/2015-05-27/land02/scripts/ggland
ggland_script=xb-sst.sh

mergeroot_dir=/n/home/bstefan/m/deexcite/git/11-05-2015/scripts/xb_runs


#quasi elastic scattering

T=7123

$egdir/quasi --target-N 0 --target-Z 1 --proj-Z 7 --proj-N 10 --cluster-Z 1 --cluster-N 0 -T $T -E 10 -J 0 -N $events> $datadir/$quasi_file
$egdir/quasi --target-N 0 --target-Z 1 --proj-Z 7 --proj-N 10 --cluster-Z 1 --cluster-N 0 -T $T -E 10.5 -J 0 -N $events>> $datadir/$quasi_file

#we do ugly things since bash cant deal very well with reals
for E in {11..20..1}
do
    echo "E: $E"
    E2=$(echo "$E +0.5" | bc)
    $egdir/quasi --target-N 0 --target-Z 1 --proj-Z 10 --proj-N 10 --cluster-Z 1 --cluster-N 0 -T $T -E $E -J 0 -N $events>> $datadir/$quasi_file
     $egdir/quasi --target-N 0 --target-Z 1 --proj-Z 10 --proj-N 10 --cluster-Z 1 --cluster-N 0 -T $T -E $E2 -J 0 -N $events>> $datadir/$quasi_file 
done


#deexcite the prefragments from the QE scattering and create a gunfile of output
$egdir/deexcite -J 10 -l 5 <$datadir/$quasi_file > $datadir/$output
$egdir/out2gun <$datadir/$output > $datadir/$gun_file
$egdir/out2root <$datadir/$output -f $datadir/$fire_file
#cp gunfile and run ggland on it
cp $datadir/$gun_file $GG
cp $datadir/$gun_file $GG
cd $GG
sleep 1
./$ggland_script $gun_file $root_file $gg_events
cp $GG/$root_file $datadir
sleep 1
cd $datadir
#add addback energies, gamma multiplicities, etc to events in the tree
./dress_tree.o