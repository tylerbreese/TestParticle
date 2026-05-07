#!/bin/bash

max_jobs=4

echo "Running bubble.exe"
echo "How many particles?"
read user_num
result=$(( (user_num + (1000 / 2)) / 1000 ))
echo "Use $result processes"

for ((i=1; i<=result; i++))
do

    while (($(jobs -r | wc -l) >= max_jobs )); do
        wait -n
    done
    echo "Running instance $i:"
    # mpirun -np 3 ./shock.exe
    ./bubble.exe & 
    sleep 30

done
wait
echo "Mission complete boss!"