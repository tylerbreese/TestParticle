./batch_multi.sh multi/input1.in & 
PID_A=$!
./batch_multi.sh multi/input2.in & 
PID_B=$!
wait  $PID_A $PID_B

./batch_multi.sh multi/input3.in &
PID_C=$! 
./batch_multi.sh multi/input4.in &
PID_D=$!
wait  $PID_C $PID_D 

./batch_multi.sh multi/input5.in &
PID_E=$!
./batch_multi.sh multi/input6.in &
PID_F=$!
wait  $PID_E $PID_F 

./batch_multi.sh multi/input7.in &
PID_G=$! 
./batch_multi.sh multi/input8.in &
PID_H=$!
wait $PID_G $PID_H 

./batch_multi.sh multi/input9.in &
PID_I=$! 
./batch_multi.sh multi/input10.in &
PID_J=$!
wait $PID_I $PID_J 