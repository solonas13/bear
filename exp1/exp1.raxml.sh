#!/bin/bash

for i in {1..50}
do
 ./raxmlHPC-PTHREADS-SSE3 -m GTRCAT -n "$i" -p 12345 -s out."$i".bear
 rm RAxML_i*
 rm RAxML_l*
 rm RAxML_p*
 rm RAxML_r*
done
