#!/bin/bash

openmpirun -np 1 ./example-opt -pc_type bjacobi
#openmpirun -np 4 ./example-opt -pc_type bjacobi
