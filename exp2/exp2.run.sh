#!/bin/bash

for i in {1..50}
do
 ../bear -a DNA -p 12.250.rot."$i".fas -k 15 -d 2 -w 15 -o out."$i".bear
done
