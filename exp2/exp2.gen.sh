#!/bin/bash

for i in {1..50}
do
 ../bear -a DNA -p 12.250.fas -k 0 -d 0 -o out -r 12.250.rot."$i".fas
done
