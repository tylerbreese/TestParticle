./batch_helium.sh helium/input1.in & 
PID_A=$!
./batch_helium.sh helium/input2.in & 
PID_B=$!
./batch_helium.sh helium/input3.in &
PID_C=$! 
./batch_helium.sh helium/input4.in &
PID_D=$!
./batch_helium.sh helium/input5.in &
PID_E=$!
wait  $PID_A $PID_B $PID_C $PID_D $PID_E 

# ./batch_helium.sh helium/input6.in &
# PID_F=$!
# ./batch_helium.sh helium/input7.in &
# PID_G=$! 
# ./batch_helium.sh helium/input8.in &
# PID_H=$!
# ./batch_helium.sh helium/input9.in &
# PID_I=$! 
# ./batch_helium.sh helium/input10.in &
# PID_J=$!
# wait $PID_F $PID_G $PID_H $PID_I $PID_J 