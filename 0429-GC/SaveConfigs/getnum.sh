#!/bin/bash
touch result.txt
for ((i=0;i<420;i++))
do 
grep 'This energy is' < RelaMin_$i.xyz | awk '{print $9 "\t" $4}' >>result.txt

done
