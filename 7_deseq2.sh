#!/bin/bash

## To generate list: for x in *gtf; do lib="${x%.gtf}"; lib="${lib#HI*In*.}"; a=$(readlink -f $x); echo -e $lib"\t"$a >> list.txt; done


prepDE.py \
    -i list.txt

