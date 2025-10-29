./batch_proton.sh proton/input1.in & 
PID_A=$!
./batch_proton.sh proton/input2.in & 
PID_B=$!
./batch_proton.sh proton/input3.in &
PID_C=$! 
./batch_proton.sh proton/input4.in &
PID_D=$!
wait  $PID_A $PID_B $PID_C $PID_D

./batch_proton.sh proton/input5.in &
PID_E=$!
./batch_proton.sh proton/input6.in &
PID_F=$!
./batch_proton.sh proton/input7.in &
PID_G=$! 
./batch_proton.sh proton/input8.in &
PID_H=$!
wait  $PID_E $PID_F $PID_G $PID_H

./batch_proton.sh proton/input9.in &
PID_I=$! 
./batch_proton.sh proton/input10.in &
PID_J=$!
wait $PID_I $PID_J