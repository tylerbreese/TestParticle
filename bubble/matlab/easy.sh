./batch_plane.sh plane/input1.in & 
PID_A=$!
./batch_plane.sh plane/input2.in & 
PID_B=$!
./batch_plane.sh plane/input3.in &
PID_C=$! 
wait  $PID_A $PID_B $PID_C 

./batch_plane.sh plane/input4.in &
PID_D=$!
./batch_plane.sh plane/input5.in &
PID_E=$!
./batch_plane.sh plane/input6.in &
PID_F=$!
wait  $PID_D $PID_E $PID_F 

./batch_plane.sh plane/input7.in &
PID_G=$! 
./batch_plane.sh plane/input8.in &
PID_H=$!
wait $PID_G $PID_H 

./batch_plane.sh plane/input9.in &
PID_I=$! 
./batch_plane.sh plane/input10.in &
PID_J=$!
wait $PID_I $PID_J 