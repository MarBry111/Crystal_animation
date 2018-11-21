#!/bin/bash
for (( i=0 ; $i <= 50 ; i++ )) ; do
python3 main.py -i $i
echo $i
done
