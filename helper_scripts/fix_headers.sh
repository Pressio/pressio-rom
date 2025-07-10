#!/usr/bin/env bash

# make sure you have guardonce installed: https://github.com/cgmb/guardonce

# this script must be run where it is
# no args needed

# rel dir
PCK_REL_DIR=${PWD}/../include/pressio

# array of names
declare -a pcks=("solvers_linear" "solvers_nonlinear" "ode" "rom")

# loop over and fix he
for packName in ${pcks[@]}; do
    # target dir
    PCK_DIR=${PCK_REL_DIR}/${packName}
    echo ${PCK_DIR}
    echo ${packName}

    cd ${PCK_DIR}
    echo ${PCK_DIR}

    # first, convert all guards to pragmas
    guard2once -r .
    # then converts from pragmas to header with specific pattern
    once2guard -r -p "path -1 | prepend PRESSIOROM_${packName}_ | upper | append _" -s "#endif  // %\n" .
done
