#!/bin/bash
touch result.txt
for ((i=0;i<610;i++))
do 
grep 'This energy is' < RelaMin_$i.xyz | awk '{print $9 "\t" $4}' >>result.txt
done
touch result_min.txt
for ((i=0;i<15;i++))
do
grep 'This energy is' < Min_$i.xyz | awk '{print $9 "\t" $4}' >>result_min.txt
done

