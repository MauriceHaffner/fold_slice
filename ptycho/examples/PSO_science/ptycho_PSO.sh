#!/bin/bash

### specify the queue: 
PBS -q defq

date

### cd to directory where the job was submitted:
cd /users/stud/haffnerm/Ptychoshelves/fold_slice/ptycho/examples/PSO_science

echo "----------------"
echo "PBS job running on: haffnerm"
echo "in directory: /users/stud/haffnerm/Ptychoshelves/fold_slice/ptycho/examples/PSO_science"
echo "----------------"

### run the program:
matlab -nodisplay -r ptycho_electron_PSO_science

date