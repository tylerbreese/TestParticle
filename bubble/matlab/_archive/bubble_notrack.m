% beta acceleration across collisionless shock boundary of bubble
% inputs for plasma speed, compression ratio, particle Z and A, and shock normal angles
% simulate performed in the shock frame using CGS units
% correlation length of turbulence fixed at k*R_E earth radii, maybe 1000km
% to start
%% inputs and constants
sample_size = 1000;
th = deg2rad(90.0); %shock angle
del = deg2rad(0.0); %streaming angle, DEL SHOULD BE FIXED AT 0, HEAD ON SHOCK
V_deb = 1400E5 .* [cos(del); 0.0; sin(del)]; %debris initial speed, cm/s
B0 = 0.5 .* [cos(th); 0.0; sin(th)]; %magnetic field in Gauss

global c e Z A q m
c = 3E10; %speed of light, m/s
e = 4.8E-10; %elementary charge, cgs
Z = 1; %atomic number
A = 1; %mass number
q = -1*e; % charge 
m = 9.11e-28; %g, mass


%% initialize
[p0,r,a,b] = sampling(sample_size,th,del,V_deb);
p = zeros(sample_size,3,3); %Velocity [Particle,Dimension(x,y,z),Time]
p(:,:,1) = p0';
p(:,:,3) = p(:,:,1);

[w_g,r_g,gamma] = calc_gyrofrequency(B0,p0');

x0 = zeros(sample_size,3); %this defines the shock location. changing moves the boundary
x1 = x0(:,1);
x1(1:sample_size/2-1) = -20.*r_g(1:sample_size/2-1);
x1(sample_size/2:end) =  20.*r_g(sample_size/2:end);
y0 = randi([-100,100],sample_size,1).*r_g;
z0 = randi([-100,100],sample_size,1).*r_g;
X = zeros(sample_size,3,2); %Position [Particle,Dimension(x,y,z),Time]
X(:,:,1) = [x0(:,1)+x1,y0,z0]; % start slighlt upstream of shock boundary
X(:,:,3) = X(:,:,1); %save initial conditions

N = 400; s = sqrt(0.5) .* vecnorm(B0,2,1); %variance s^2 / B^2 = number < 1
[An,kk] = init_turby(V_deb,w_g,s,N); %turby on
%An = zeros(sample_size,N); kk=An; %turby off
[B,U] = shock_field(B0,X(:,:,1),r,a,b,x0,V_deb,sample_size,An,kk);

%% integration
dt = 0.05 * (1/w_g(1)); %might need smaller timestep for shock frame
n = 1; t = 0.0;
Split = {};
%while t < 0.0034 %30,000 cyclotron periods
while n < 300000*20+1 %300k cyclotron periods
% if n > 2000
%     break
% end

    % %BORIS
    % E(:,1) = (-1/c) .* (U(:,2).*B(:,3) - U(:,3).*B(:,2));
    % E(:,2) = (-1/c) .* (U(:,3).*B(:,1) - U(:,1).*B(:,3));
    % E(:,3) = (-1/c) .* (U(:,1).*B(:,2) - U(:,2).*B(:,1));
    % 
    % R = (q) .* E .* (dt./2); 
    % 
    % Pm = p(:,:,1) + R;
    % %gamma = sqrt( 1 + sum(Pm.*Pm,2)./(m^2 * c^2) );
    % T = (q./(gamma.*(m*c))) .* B .* (dt./2);
    % S = 2.*T ./ (1 + sum(T.*T,2));
    % 
    % Pp(:,1) = Pm(:,1) + ( Pm(:,2).*T(:,3) - Pm(:,3).*T(:,2) );
    % Pp(:,2) = Pm(:,2) + ( Pm(:,3).*T(:,1) - Pm(:,1).*T(:,3) );
    % Pp(:,3) = Pm(:,3) + ( Pm(:,1).*T(:,2) - Pm(:,2).*T(:,1) );
    % 
    % Pq(:,1) = Pm(:,1) + ( Pp(:,2).*S(:,3) - Pp(:,3).*S(:,2) );
    % Pq(:,2) = Pm(:,2) + ( Pp(:,3).*S(:,1) - Pp(:,1).*S(:,3) );
    % Pq(:,3) = Pm(:,3) + ( Pp(:,1).*S(:,2) - Pp(:,2).*S(:,1) );
    % 
    % p(:,:,2) = Pq + R;
    % 
    % X(:,1,2) = X(:,1,1) + (p(:,1,2)./(gamma.*m)) .* dt;
    % X(:,2,2) = X(:,2,1) + (p(:,2,2)./(gamma.*m)) .* dt;
    % X(:,3,2) = X(:,3,1) + (p(:,3,2)./(gamma.*m)) .* dt;
    %integrate
    [X(:,:,2),p(:,:,2)] = integrate_boris(X(:,:,1),p(:,:,1),U,B,dt,gamma);  
    cutoff = logspace(log10(2),1,10);
    [Split] = particle_split(X,p,cutoff,Split);
    [Split] = particle_displace(Split,vecnorm(V_deb,2,1)/w_g);
    for ii = 1:length(Split)
        if isempty(Split{ii})
            [X_split,V_split] = integrate_boris(Split{ii}{1},Split{ii}{2},U,B,dt);
            Split{ii}{1} = X_split; Split{ii}{2} = V_split;
        else
            continue
        end
    end
    n = n + 1; t = t + dt; 
    %advance field
    [B,U] = shock_field(B0,X(:,:,2),r,a,b,x0,V_deb,sample_size,An,kk);
 
    [V_deb]  = get_debris_speed(t,dt,del);
    %[r,a,b] = get_compression_ratio(V_deb,th,del);
    [w_g,r_g,gamma] = calc_gyrofrequency(B0,p(:,:,2));
    %scattering
    %[X,p(:,:,2)] = scatter_atmos(X,p(:,:,2),gamma);
    %
    X(:,:,1) = X(:,:,2);
    p(:,:,1) = p(:,:,2); 
    n = n + 1; t = t + dt; 
end
%sim end save data
h5create('output_shock.h5','/position',size(X));
h5write('output_shock.h5','/position',X);
h5create('output_shock.h5','/momentum',size(p));
h5write('output_shock.h5','/momentum',p);
h5create('output_shock.h5','/gyro',[1,2]);
h5write('output_shock.h5','/gyro',[r_g, w_g.*ones(sample_size,1)]);
% h5create('output.h5','/particle',[1,4]);
% h5write('output.h5','/particle',[m,q,Z,A]);

%% functions
function [B,U] = shock_field(B0,X,r,a,b,x0,U0,sample_size,An,kk)
    % for the nth time cycle update B based on location 
    %if upstream of shock B0
    %if downstream of shock r*B0
    %shock taken in xz plane with downstream being the positive direction
    track = double(X(:,1) > x0(:,1));
    
    U = ones(sample_size,3);
 
    U(track==0,1) = U0(1);
    U(track==0,2) = U0(2);
    U(track==0,3) = U0(3);
    U(track~=0,1) = U0(1)./r;
    U(track~=0,2) = U0(2) + b .* U0(1);
    U(track~=0,3) = U0(3) + b .* U0(1);
    
    B = ones(sample_size,3);
    N = width(kk); [dB,~] = sum_turby(X,N,kk,An);
    %dB = zeros(size(B)); % turby off
    B(track==0,1) = B0(1) + dB(track==0,1);
    B(track==0,2) = B0(2) + dB(track==0,2);
    B(track==0,3) = B0(3) + dB(track==0,3);
    B(track~=0,1) = B0(1) + dB(track~=0,1);
    B(track~=0,2) = ( B0(2) + dB(track~=0,2) );
    B(track~=0,3) = a .* ( B0(3) + dB(track~=0,3) );
        
end

function [X] = pdfrnd(x, px, sampleSize)
    cdf = cumsum(px);
    cdf = cdf/sum(cdf);
    cdf = cdf/max(cdf);
    rnd = rand(sampleSize, 1);
    X = interp1(cdf, x, rnd, 'linear', 0);
end

function [trapped,escaped] = check_escape(Vmag,V_deb,r)
U = norm(V_deb)./r;

%P = ones(size(Vmag)) - ( ((Vmag - U).^2) ./ ((Vmag + U).^2) );
%check = rand(size(Vmag));
%track = double(check < P); %two different methods
track = double(Vmag<U); %this gives more particles

trapped = zeros(size(Vmag));
trapped(track==0) = 1;
trapped(track~=0) = 0;

escaped = zeros(size(Vmag));
escaped(track~=0) = 1;
escaped(track==0) = 0;

end

function [An,kk] = init_turby(U,w_g,s,N)
    %initialize spectrum
    %needs position, local streaming velocity and frequency for each trajectory
    %user pick Number of modes
    r_gn = abs(vecnorm(U,2,1) ./ w_g(1));
    L = 10*1e5; %1000km in cm
    l1 = 5 * L; l2 = 0.5*r_gn;
    [kk,dkk] = get_k_vec(l2,l1,N);
    
    gamma = 11/3; %3D turb`
    P0 = kk.^2 .* dkk; %make the psd array
    GkN = P0 ./ (1 + (kk.*L).^gamma);
    Gk = sum(GkN,2);
    
    An = 2.0 * sqrt( s^2 ./ Gk) .* sqrt(GkN); 
    
end
    
function [K,dK] = get_k_vec(k_min,k_max,N)
    %k goes as U/w_g
    %returns n by N matrix of N k vectors for n samples
    %logspace between kmin and kmax
    % excessive but i like this in a subroutine
    %N+1 so you get same modes as diff
    n = length(k_max);
    K = zeros(n,N+1);
        for ii = 1:n
            a = log10(k_min(ii)); b = log10(k_max(ii));
            K(ii,:) = logspace(a,b,N+1);
        end
    %get delta K
    dK = diff(K,1,2); K = K(:,2:end);
end

function [w_g,r_g,gamma] = calc_gyrofrequency(B,p)
    global c q m 
    gamma = (1 + ((vecnorm(p,2,2)./(m*c))).^2).^(0.5);
   
    % 
    % boost(:,1) = (1 - (p(:,1).^2 ./ m*c^2))./gamma;
    % boost(:,2) = (1 - (p(:,2).^2 ./ m*c^2))./gamma;
    % boost(:,3) = (1 - (p(:,3).^2 ./ m*c^2))./gamma;
    % 
    %w_g = (q/m) .* (vecnorm(B,2,1)./gamma);
    w_g = (q .* vecnorm(B,2,1)) ./  ( m * c);
    w_g = abs(w_g);
    r_g = vecnorm(p,2,2) ./ (q * vecnorm(B,2,1));
    r_g = abs(r_g);
end

function [dB,dBn] = sum_turby(X,N,kk,An)

    persistent turby

    if isempty(turby)
        [n,~,~] = size(X); %grab sample size
        beta = (2*pi).*rand(N,1); %random phase
        u = 2.*rand(n,1) - 1; % rand on -1 to 1
        %th = acos(u);
        phi = (2*pi).*rand(n,1);
        alpha = (2*pi).*rand(n,1);
        
        % Calculate rotation matrix components
        sin_u = sqrt(1 - u.^2);  % sqrt(1 - u^2)
        
        %make wavelets
        Anx =  An .* cos(alpha) .* cos(phi) .* u;
        Any =  An .* cos(alpha) .* sin(phi) .* u;
        Anz = -An .* cos(alpha) .* sin_u;
        Bnx =  An .* sin(alpha) .* sin(phi);
        Bny = -An .* sin(alpha) .* cos(phi);
    
        kx = kk .* cos(phi) .* sin_u;
        ky = kk .* sin(phi) .* sin_u;
        kz = kk .* u;

    % Store everything in persistent struct
    turby.kx = kx; turby.ky = ky; turby.kz = kz;
    turby.beta = beta;
    turby.Anx = Anx; turby.Bnx = Bnx;
    turby.Any = Any; turby.Bny = Bny;
    turby.Anz = Anz;
    end

    arg = (turby.kx .* X(:,1) + turby.ky .* X(:,2)+ turby.kz .* X(:,3) + turby.beta');
    % sum
    dBn(:,:,1) = turby.Anx .* cos(arg) + turby.Bnx .* sin(arg);
    dBn(:,:,2) = turby.Any .* cos(arg) + turby.Bny .* sin(arg);
    dBn(:,:,3) = turby.Anz .* cos(arg);
    dB = squeeze(sum(dBn,2));
end


%% sampling
function [p0,r,a,b] = sampling(sample_size,th,del,V_deb)
    global m c
    [r,a,b] = get_compression_ratio(V_deb,th,del);
    %beta energy dist, assume U235
    in = readmatrix('THTk_default_betas.txt');
    CR = in(3:203,1:2);
    E_plasma = pdfrnd(CR(:,1),CR(:,2),sample_size);
    E_plasma = 1.6e-6 .* E_plasma; %convert to ergs
    g = sqrt(1 + (E_plasma) ./ (m*c*c));
    p_plasma = (m*c) .* (g.^2 -1).^(0.5);
    %spherical particle distribution
    u = pdfrnd(-1:0.001:1,ones(size(-1:0.001:1)),sample_size);
    phi_max = 2*pi; phi_min = 0;
    div = abs(phi_max-phi_min)/length(u); %force square array of samples
    phi = pdfrnd(phi_min:div:phi_max,ones(size(phi_min:div:phi_max)),sample_size);
    p0 = [
        (p_plasma .* ( 1 - u.^2).^0.5 .* cos(phi))'
        (p_plasma .* ( 1 - u.^2).^0.5 .* sin(phi))'
        (p_plasma .* u)'
        ];
    % move into shock frame maybe doesn't matter
    %p0 = p0 + m.*V_deb; %
    
    %try check to replace vectors where V_x0 <0. onlky want forward moving particles; 
    % for ii = 1:sample_size
    %     while V0(1,ii) < 0
    %         u = randsample(seed,-1:0.001:1,1);
    %         phi_max = 2*pi; phi_min = 0;
    %         div = abs(phi_max-phi_min)/length(u); %force square array of samples
    %         phi = randsample(seed,phi_min:div:phi_max,1);
    %         V0(:,ii) =  [
    %             V_plasma(ii) .* ( 1 - u.^2).^0.5 .* cos(phi)
    %             V_plasma(ii) .* ( 1 - u.^2).^0.5 .* sin(phi)
    %             V_plasma(ii) .* u
    %             ]; %resample until V0_x is moving forward
    %     end
    % end
end

function [U] = get_debris_speed(t,dt,del)

    U0 = 1400E5 .* [cos(del); 0.0; sin(del)]; %just leave this fixed
    % Assign the fixed debris speed to U
    U = U0; %constant (colgate 0th order)
    %U = U0 .* (t/dt).^(-1); %stupisky
    %U = U0 .* (t/dt).^(-3/5); % taylor sedov idk maybe fun to try

end
function [r,a,b] = get_compression_ratio(U,th,del)
    Va = 40e5; %cm/s
    Ma = vecnorm(U,2,1) ./ Va;
    r = 6; %from dyal
    a = r * (Ma^2*cos(del)^2 - cos(th)^2) / (Ma^2*cos(del)^2-r*cos(th)^2);
    b = ((r-1)*cos(th)*sin(th)) ./ (Ma^2*cos(del)^2-r*cos(th)^2);
end

function [X_scat,p_scat] = scatter_atmos(Xn,pn,gamma)

Zeff = 8.0; %atomic oxygen dominate species
a0 = 5.29e-9; %bohr radius
a1 = 0.885*a0 / (Zeff)^(1/3);

K = 4*pi * ((1/137)*Zeff*a1)^2;

beta = sqrt(1 - 1./gamma.^2);

sigma = K ./ beta.^2;
nd = 16*6.022e23* 1e-16;
L = 1./(nd.*sigma);

ds = sqrt(Xn(:,1,2).^2 + Xn(:,2,2).^2 + Xn(:,3,2).^2) - sqrt(Xn(:,1,1).^2 + Xn(:,2,1).^2 + Xn(:,3,1).^2);

P    = 1 - (abs(ds)./L);
roll = rand(size(ds));

track = double(roll < P);
%track = rand(100, 1) < 0.25; %force
n = length(nonzeros(track));
u = 2.*rand(n,1) - 1; % rand on -1 to 1
%th = acos(u);
phi = (2*pi).*rand(n,1);
% Calculate rotation matrix components
sin_phi = sin(phi);
cos_phi = cos(phi);
sin_u = sqrt(1 - u.^2);  % sqrt(1 - u^2)

p_scat = ones(size(pn));
p_scat(track==0,1) = pn(track==0,1);
p_scat(track==0,2) = pn(track==0,2);
p_scat(track==0,3) = pn(track==0,3);

R = zeros(n, 3, 3);  % Initialize the 3D matrix to store the rotation matrices

R(:, 1, 1) = u .* cos_phi;  % First column of R first row
R(:, 1, 2) = u .* sin_phi;  % Second column of R
R(:, 1, 3) = sin_u;        % Third column of R

R(:, 2, 1) = -sin_phi;  % First column of second row of R
R(:, 2, 2) = cos_phi;   % Second column of second row of R
R(:, 2, 3) = 0;         % Third column of second row of R

R(:, 3, 1) = -sin_u .* cos_phi;  % First column of third row of R
R(:, 3, 2) = -sin_u .* sin_phi;  % Second column of third row of R
R(:, 3, 3) = u;             % Third column of third row of R

p_scat(track==1,:) = reshape(sum(R .* reshape(pn(track==1,:), [n, 1, 3]), 3), n, 3);
% p_scat(track==1,1) = pn(track==1,1) .* sin_u .* cos_phi;
% p_scat(track==1,2) = pn(track==1,2) .* sin_u .* sin_phi;
% p_scat(track==1,3) = pn(track==1,3) .* u;

X_scat = Xn;

end


function [X,p] = integrate_boris(Xn,pn,U,B,dt,gamma)
    global q m c

    %BORIS
    E(:,1) = (-1/c) .* (U(:,2).*B(:,3) - U(:,3).*B(:,2));
    E(:,2) = (-1/c) .* (U(:,3).*B(:,1) - U(:,1).*B(:,3));
    E(:,3) = (-1/c) .* (U(:,1).*B(:,2) - U(:,2).*B(:,1));

    R = (q) .* E .* (dt./2); 
    
    Pm = pn + R;
    %gamma = sqrt( 1 + sum(Pm.*Pm,2)./(m^2 * c^2) );
    T = (q./(gamma.*(m*c))) .* B .* (dt./2);
    S = 2.*T ./ (1 + sum(T.*T,2));

    Pp(:,1) = Pm(:,1) + ( Pm(:,2).*T(:,3) - Pm(:,3).*T(:,2) );
    Pp(:,2) = Pm(:,2) + ( Pm(:,3).*T(:,1) - Pm(:,1).*T(:,3) );
    Pp(:,3) = Pm(:,3) + ( Pm(:,1).*T(:,2) - Pm(:,2).*T(:,1) );

    Pq(:,1) = Pm(:,1) + ( Pp(:,2).*S(:,3) - Pp(:,3).*S(:,2) );
    Pq(:,2) = Pm(:,2) + ( Pp(:,3).*S(:,1) - Pp(:,1).*S(:,3) );
    Pq(:,3) = Pm(:,3) + ( Pp(:,1).*S(:,2) - Pp(:,2).*S(:,1) );

    p       = Pq + R;
    
    X(:,1) = Xn(:,1) + (p(:,1)./(gamma.*m)) .* dt;
    X(:,2) = Xn(:,2) + (p(:,2)./(gamma.*m)) .* dt;
    X(:,3) = Xn(:,3) + (p(:,3)./(gamma.*m)) .* dt;
end  

function [split] = particle_split(X,p,cutoff,split)
    global m c
    En = sqrt(vecnorm(p(:,:,2),2,2).^2 .* c^2 + (m*c^2)^2);
    E0 = sqrt(vecnorm(p(:,:,3),2,2).^2 .* c^2 + (m*c^2)^2);
    ratio = En ./ E0;
    flag = ratio>cutoff;
    if isempty(split)
        for ii = 1:length(cutoff)
            split{ii} = {X(flag(:,ii),:,2), p(flag(:,ii),:,2), flag(:,ii)};
        end
    else
        temp_flag = flag;
        for ii = 1:length(split)
            flag(:,ii) = split{ii}{3};
        end
        flag = flag | temp_flag;
        flag(flag>1) = 1;
        for ii = 1:length(cutoff)
            ssplit{ii} = {X(flag(:,ii),:,2), p(flag(:,ii),:,2), flag(:,ii)};
        end
        split = ssplit;
    end    
end

function [split] = particle_displace(split,R_g)
    if ~isempty(split)
        for ii = 1:length(split)
            X_split = split{ii}{1};
            X_split(:,1) = X_split(:,1) + rand(height(X_split),1).*R_g;
            X_split(:,2) = X_split(:,2) + rand(height(X_split),1).*R_g;
            X_split(:,3) = X_split(:,3) + rand(height(X_split),1).*R_g;
            V_split = split{ii}{2};
            V_split = rand_rotate(V_split);
            split{ii}{1} = X_split; split{ii}{2} = V_split;
        end
    end
end


function [Xp] = rand_rotate(X)
    [n,~,~] = size(X); %grab sample size
    u = 2.*rand(n,1) - 1; % rand on -1 to 1
    phi = (2*pi).*rand(n,1);
    % Calculate rotation matrix components
    sin_phi = sin(phi);
    cos_phi = cos(phi);
    sin_u = sqrt(1 - u.^2);  % sqrt(1 - u^2)
    
    % Build the rotation matrix R for all rows at once
    R = zeros(n, 3, 3);  % Initialize the 3D matrix to store the rotation matrices
    
    R(:, 1, 1) = u .* cos_phi;  % First column of R first row
    R(:, 1, 2) = u .* sin_phi;  % Second column of R
    R(:, 1, 3) = sin_u;        % Third column of R
    
    R(:, 2, 1) = -sin_phi;  % First column of second row of R
    R(:, 2, 2) = cos_phi;   % Second column of second row of R
    R(:, 2, 3) = 0;         % Third column of second row of R
    
    R(:, 3, 1) = -sin_u .* cos_phi;  % First column of third row of R
    R(:, 3, 2) = -sin_u .* sin_phi;  % Second column of third row of R
    R(:, 3, 3) = u;             % Third column of third row of R
    
    % Rotate X
    Xp = reshape(sum(R .* reshape(X, [n, 1, 3]), 3), n, 3);
end
