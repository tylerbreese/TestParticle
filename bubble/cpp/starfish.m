% calc plasma parameters and dispersion relations for starfish 

%% fields
B  = 0.5; % geomagnetic field
U  = 2000e5; %shock speed
vA = 300e5; %alfven speed
Cs = 0.6e5; %sound speed
T  = 1.2e7; %1 keV temp in Kelvin

y    = 7/5; %adiabatic index of air
beta = (2/y) * ( (U/vA)/(U/Cs) )^2; 