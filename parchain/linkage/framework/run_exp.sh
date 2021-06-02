#!/bin/bash

# $1: file name
# $2: 1 if want to re compile

if [ $3 -eq 1 ]
then
    make clean
    make -j
fi

./linkage -d $2 -method $4 ../completeLinkage/datasets/$1.pbbs