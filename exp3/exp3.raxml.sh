#!/bin/bash

 ./raxmlHPC-PTHREADS-SSE3 -m GTRCAT -n 50.2500.35.true -p 12345 -s 50.2500.35.aln.true
 rm RAxML_i*
 rm RAxML_l*
 rm RAxML_p*
 rm RAxML_r*
 ./raxmlHPC-PTHREADS-SSE3 -m GTRCAT -n 50.2500.35.nobear -p 12345 -s 50.2500.35.aln.nobear
 rm RAxML_i*
 rm RAxML_l*
 rm RAxML_p*
 rm RAxML_r*
 ./raxmlHPC-PTHREADS-SSE3 -m GTRCAT -n 50.2500.35.bear -p 12345 -s 50.2500.35.aln.bear
 rm RAxML_i*
 rm RAxML_l*
 rm RAxML_p*
 rm RAxML_r*

 cat RAxML_bestTree.50.2500.35.true RAxML_bestTree.50.2500.35.bear RAxML_bestTree.50.2500.35.nobear > 50.2500.35.allTrees
 ./raxmlHPC-PTHREADS-SSE3 -f r -z 50.2500.35.allTrees -m GTRCAT -n 50.2500.35
 rm RAxML_i*
