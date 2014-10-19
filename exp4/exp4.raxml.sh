#!/bin/bash

 ./raxmlHPC-PTHREADS-SSE3 -m GTRCAT -n primates.bear -p 12345 -s primates.bear.aln
 rm RAxML_i*
 rm RAxML_l*
 rm RAxML_p*
 rm RAxML_r*
 ./raxmlHPC-PTHREADS-SSE3 -m GTRCAT -n primates.cyc -p 12345 -s primates.cyc.aln
 rm RAxML_i*
 rm RAxML_l*
 rm RAxML_p*
 rm RAxML_r*

 cat RAxML_bestTree.primates.cyc RAxML_bestTree.primates.bear > primates.allTrees
 ./raxmlHPC-PTHREADS-SSE3 -f r -z primates.allTrees -m GTRCAT -n primates.allTrees
 rm RAxML_i*
