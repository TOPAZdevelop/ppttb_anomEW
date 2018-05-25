#!/bin/bash


# generate anomalous coupling histograms

List1='-05.00 -02.00 +00.00 +02.00 +05.00'

for Cphiq in $List1; do  
for Cphiu in $List1; do  


    echo "($Cphiq,$Cphiu)"; 
    
    ./ReadCSV ./input_LO.dat ./input_couplings.dat $Cphiq $Cphiu CouplingData/ttbarEW_${Cphiq}_${Cphiu}.dat

done
done

 
 

 
# run LL analysis program


 for Cphiq in $List1; do  
for Cphiu in $List1; do  


    echo "($Cphiq,$Cphiu)"; 
    mpirun -n 24 ../../BinnedLogL2 GetLogL ./CouplingData/ttbarEW_+00.00_+00.00.dat CouplingData/ttbarEW_${Cphiq}_${Cphiu}.dat 1 300 10000 5.0 AnalysisData/ttbarEW_${Cphiq}_${Cphiu}.dat

done
done
