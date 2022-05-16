#!/bin/bash 


Dir="/groups/pupko/cohenofi/HGT/runCOG.iTOL.38.3915/Batch.V9.93.PostExpSim.FewMPcostRatios.test"
c="4"
s=38
g=3915


for i in 1 2 3 4 5 6 7 8 9 10
do


for d in  GL$c  GL$c.MP1params GL$c.MP2params GL$c.MP3params GL$c.FixedParams GL$c.FixedParams
do

for r in 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5
do

echo "Perl readPostExpSimulation $Dir/AcountForMissingData.$d/RESULTS/SimulatedPostExp$i  MP with cost $r"
perl /groups/pupko/cohenofi/pupkoSVN/trunk/scripts/gainLoss/readPostExpSimulation.pl  $Dir/AcountForMissingData.$d/RESULTS/SimulatedPostExp$i/ T $s $g $r
cp  $Dir/AcountForMissingData.$d/RESULTS/SimulatedPostExp$i/outRes.txt  $Dir/AcountForMissingData.$d/RESULTS/
mv $Dir/AcountForMissingData.$d/RESULTS/outRes.txt  $Dir/AcountForMissingData.$d/RESULTS/outResMP$i.$r.txt

done


echo "Perl readPostExpSimulation $Dir/AcountForMissingData.$d/RESULTS/SimulatedPostExp$i"
perl /groups/pupko/cohenofi/pupkoSVN/trunk/scripts/gainLoss/readPostExpSimulation.pl  $Dir/AcountForMissingData.$d/RESULTS/SimulatedPostExp$i/
cp  $Dir/AcountForMissingData.$d/RESULTS/SimulatedPostExp$i/outRes.txt  $Dir/AcountForMissingData.$d/RESULTS/
mv $Dir/AcountForMissingData.$d/RESULTS/outRes.txt  $Dir/AcountForMissingData.$d/RESULTS/outRes$i.txt

echo "Perl readPostExpSimulationPerPos $Dir/AcountForMissingData.$d/RESULTS/SimulatedPostExp$i"
perl /groups/pupko/cohenofi/pupkoSVN/trunk/scripts/gainLoss/readPostExpSimulationPerPos.pl $Dir/AcountForMissingData.$d/RESULTS/SimulatedPostExp$i/
cp  $Dir/AcountForMissingData.$d/RESULTS/SimulatedPostExp$i/outResPerPos.txt  $Dir/AcountForMissingData.$d/RESULTS/
mv $Dir/AcountForMissingData.$d/RESULTS/outResPerPos.txt  $Dir/AcountForMissingData.$d/RESULTS/outResPerPos$i.txt

echo "tar $Dir/AcountForMissingData.$d/RESULTS/SimulatedPostExp$i"
tar -cf $Dir/AcountForMissingData.$d/RESULTS/SimulatedPostExp$i.tar $Dir/AcountForMissingData.$d/RESULTS/SimulatedPostExp$i
echo "rm $Dir/AcountForMissingData.$d/RESULTS/SimulatedPostExp$i"
rm -R $Dir/AcountForMissingData.$d/RESULTS/SimulatedPostExp$i


done

done
