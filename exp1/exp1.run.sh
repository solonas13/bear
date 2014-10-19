#!/bin/bash

for i in {1..50}
do
 ../bear -a DNA -p 12.250.rot."$i".fas -k 45 -d 1 -o out."$i".bear
done
