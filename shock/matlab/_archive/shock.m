tic
% ion acceleration across collisionless shock boundary
% inputs for plasma speed, compression ratio, particle Z and A, and shock normal angles
% simulate performed in the shock frame using CGS units
% correlation length of turbulence fixed at 1AU
%% inputs and constants
sample_size = 200;
%dt = 1; %secs ~0.05 (1/w_g)
cycles = 300; %number of time segments
r = 2.5; %compression ratio, %strong r>4, weak shock ~1.25<=r<4
th = deg2rad(90.0); %shock angle
del = deg2rad(0.0); %streaming angle, DEL SHOULD BE FIXED AT 0, HEAD ON SHOCK
V_SW = 400E5 .* [cos(del); 0.0; sin(del)]; %solar wind speed, cm/s
B0 = 0.05E-5 .* [cos(th); 0.0; sin(th)]; %magnetic field in Tesla B(r) = B0/R[AU]
% %read in file
% fid         = fopen('run_072325_VSW_400_r_4.in','r');
% formatSpec  = '%s %s %s\n';
% inputdata   = textscan(fid,formatSpec,8);
% run_name    = join([string(inputdata{1,3}{1}),'.h5'],'None');
% sample_size = str2double(inputdata{1,3}{2});
% cycles      = str2double(inputdata{1,3}{3});
% r           = str2double(inputdata{1,3}{4});
% th          = str2double(inputdata{1,3}{5});
% del         = str2double(inputdata{1,3}{6});
% VSW         = str2double(inputdata{1,3}{7});
% B_0         = str2double(inputdata{1,3}{8});
% V_SW        = VSW .* [cos(del); 0.0; sin(del)]; %solar wind speed, cm/s
% B0          = B_0 .* [cos(th); 0.0; sin(th)]; %magnetic field in Tesla B(r) = B0/R[AU]

global c e Z A q m
c = 3E10; %speed of light, m/s
e = 4.8E-10; %elementary charge, cgs
Z = 2; %atomic number
A = 4; %mass number
q = 1*e; %ion charge 
m = A * 1.67E-24; %g

%% initialize
[V0,a,b] = sampling(sample_size,r,th,del,V_SW,B0);
V = zeros(sample_size,3,cycles); %Velocity [Particle,Dimension(x,y,z),Time]
V(:,:,1) = V0';

[w_g,r_g,gamma,boost] = calc_gyrofrequency(B0,V0');

U1 = vecnorm(V_SW,2,1);
B1 = vecnorm(B0,2,1);
Om = q*B1 ./ (m*c);
Rg = U1 ./ Om;
x0 = zeros(sample_size,3); %this defines the shock location. changing moves the boundary
y0 = randi([-500,500],sample_size,1).*Rg;
z0 = randi([-500,500],sample_size,1).*Rg;
X = zeros(sample_size,3,cycles); %Position [Particle,Dimension(x,y,z),Time]
X(:,:,1) = [x0(:,1)-20.*Rg,y0,z0]; % start slighlt upstream of shock boundary

N = 400; s = sqrt(0.5) .* vecnorm(B0,2,1); %variance s^2 / B^2 = number < 1
[An,kk] = init_turby(V_SW,Om,s,N); %turby on
%An = zeros(sample_size,N); kk=An; %turby off
[B,U] = shock_field(B0,X(:,:,1),r,a,b,x0,V_SW,sample_size,An,kk);

Vmag(:,1) = vecnorm(V0,2,1);

%% integration
dt = 0.05 * (1/w_g(1)); %might need smaller timestep for shock frame
n = 1; t = 0.0;

while t < cycles*(1/w_g(1))
    
    E(:,1) = (-1/c) .* (U(:,2).*B(:,3) - U(:,3).*B(:,2));
    E(:,2) = (-1/c) .* (U(:,3).*B(:,1) - U(:,1).*B(:,3));
    E(:,3) = (-1/c) .* (U(:,1).*B(:,2) - U(:,2).*B(:,1));
    %BORIS
    u = V(:,:,n);
    R = (q/m) .* E .* (dt./2); 
    T = (q/(m*c)) .* B .* (dt./2);
    S = 2.*T ./ (1 + sum(T.*T,2));

    u = u + R; 

    Vp(:,1) = u(:,1) + ( u(:,2).*T(:,3) - u(:,3).*T(:,2) );
    Vp(:,2) = u(:,2) + ( u(:,3).*T(:,1) - u(:,1).*T(:,3) );
    Vp(:,3) = u(:,3) + ( u(:,1).*T(:,2) - u(:,2).*T(:,1) );

    Vq(:,1) = u(:,1) + ( Vp(:,2).*S(:,3) - Vp(:,3).*S(:,2) );
    Vq(:,2) = u(:,2) + ( Vp(:,3).*S(:,1) - Vp(:,1).*S(:,3) );
    Vq(:,3) = u(:,3) + ( Vp(:,1).*S(:,2) - Vp(:,2).*S(:,1) );

    V(:,:,n+1) = Vq + R;

    X(:,1,n+1) = X(:,1,n) + V(:,1,n+1) .* dt;
    X(:,2,n+1) = X(:,2,n) + V(:,2,n+1) .* dt;
    X(:,3,n+1) = X(:,3,n) + V(:,3,n+1) .* dt;
       
    %advance field
    [B(:,:,n),U(:,:,n)] = shock_field(B0,X(:,:,n),r,a,b,x0,V_SW,sample_size,An,kk);
    %[B,U] = shock_field(B0,X(:,:,n),r,a,b,x0,V_SW,sample_size,An,kk);
    %[~,~,gamma,boost] = calc_gyrofrequency(B',V(:,:,n+1));
    Vmag(:,n) = sqrt( V(:,1,n+1).^2 + V(:,2,n+1).^2 + V(:,3,n+1).^2 );
%    [trapped,escaped] = check_escape(Vmag,V_SW,r);

    n = n + 1; t = t + dt; 

end
%sim end save data
% h5create('output.h5','/position',size(X));
% h5write('output.h5','/position',X);
% h5create('output.h5','/velocity',size(V));
% h5write('output.h5','/velocity',V);
% h5create('output.h5','/gyro',[sample_size,2]);
% h5write('output.h5','/gyro',[r_g, w_g]);
% h5create('output.h5','/particle',[1,4]);
% h5write('output.h5','/particle',[m,q,Z,A]);
toc
run_name = 'output.h5';
h5create(run_name,'/position',size(X));
h5write(run_name,'/position',X);
h5create(run_name,'/velocity',size(V));
h5write(run_name,'/velocity',V);
h5create(run_name,'/gyro',[sample_size,2]);
h5write(run_name,'/gyro',[r_g, w_g]);
h5create(run_name,'/particle',[1,4]);
h5write(run_name,'/particle',[m,q,Z,A]);


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
    B(track~=0,2) = a .* ( B0(2) + dB(track~=0,2) );
    B(track~=0,3) = a .* ( B0(3) + dB(track~=0,3) );
        
end

function [X] = pdfrnd(x, px, sampleSize)
    cdf = cumsum(px);
    cdf = cdf/sum(cdf);
    cdf = cdf/max(cdf);
    rnd = rand(sampleSize, 1);
    X = interp1(cdf, x, rnd, 'linear', 0);
end

function [trapped,escaped] = check_escape(Vmag,V_SW,r)
U = norm(V_SW)./r;

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
    r_gn = vecnorm(U,2,1) ./ w_g(1);
    L = 1.496e+13; %1 AU in cm
    l1 = 5 * L; l2 = 0.5*r_gn;
    [kk,dkk] = get_k_vec(l2,l1,N);
    
    gamma = 11/3; %3D turb`
    P0 = kk.^2 .* dkk; %make the psd array
    GkN = P0 ./ (1 + (kk.*L).^gamma);
    Gk = sum(GkN,2);
    
    An = 2.0 * sqrt( s^2 ./ Gk) .* sqrt(GkN); 
end
    
function [K,dK] = get_k_vec(l_min,l_max,N)
    %returns n by N matrix of N k vectors for n samples
    %logspace between kmin and kmax
    % excessive but i like this in a subroutine
    %N+1 so you get same modes as diff
    k_min = (2*pi)/l_max;
    k_max = (2*pi)/l_min;
    n = length(l_max);
    K = zeros(n,N+1);
        for ii = 1:n
            a = log10(k_min(ii)); b = log10(k_max(ii));
            K(ii,:) = logspace(a,b,N+1);
        end
    %get delta K
    dK = diff(K,1,2); K = K(:,2:end);
end

function [w_g,r_g,gamma,boost] = calc_gyrofrequency(B,V)
    global c q m 
    %V = X(4:6,:) ./ m;
    Vmag = sqrt( V(:,1).^2 + V(:,2).^2 + V(:,3).^2 );
    %Vmag = vecnorm(V,2,2);
    beta = Vmag ./ c;
    gamma = 1 ./ sqrt(1 - beta.^2);
    
    boost(:,1) = (1 - (V(:,1).^2 ./ c^2))./gamma;
    boost(:,2) = (1 - (V(:,2).^2 ./ c^2))./gamma;
    boost(:,3) = (1 - (V(:,3).^2 ./ c^2))./gamma;
    
    %w_g = (q/m) .* (vecnorm(B,2,1)./gamma);
    w_g = (q .* vecnorm(B,2,1)) ./  (gamma * m * c);
    r_g = Vmag ./ w_g;
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

function [V0,a,b] = sampling(sample_size,r,th,del,V_SW,B0)
    %% sampling
    Va = 3e6; %cm/s
    Ma = vecnorm(V_SW,2,1) ./ Va;
    a = r * (Ma^2*cos(del)^2 - cos(th)^2) / (Ma^2*cos(del)^2-r*cos(th)^2);
    b = ((r-1)*cos(th)*sin(th)) ./ (Ma^2*cos(del)^2-r*cos(th)^2);
    
    %spherical particle distribution
    v = 1000:1000:vecnorm(V_SW,2,1);
    u = pdfrnd(-1:0.001:1,ones(size(-1:0.001:1)),sample_size);
    phi_max = 2*pi; phi_min = 0;
    div = abs(phi_max-phi_min)/length(u); %force square array of samples
    phi = pdfrnd(phi_min:div:phi_max,ones(size(phi_min:div:phi_max)),sample_size);
    %phi = phi'; u = u'; %transpose
    % %Siscoe
    % %C = (3/(8/pi())) * (N0*v0/(V_SW^4));
    % R = 1/100;
    % f = @(X,Y) (vecnorm(V_SW,2,1) ./ X).^(3/4) .* R .* exp(-R .* (Y./sin(Y)) .* (vecnorm(V_SW,2,1)./X).^(3/2) ); %V&siscoe dist function
    % F = f(v,pi/2);
    % V_plasma = pdfrnd(v,F,sample_size
    % V0 = [
    %   V_plasma .* ( 1 - u.^2).^0.5 .* cos(phi)
    %   V_plasma .* ( 1 - u.^2).^0.5 .* sin(phi)
    %   V_plasma .* u
    % ];
    % 
    %power law
    power1 = 1; F = v.^(power1); 
    F = F ./ trapz(v,F);
    V_plasma = pdfrnd(v,F,sample_size);
    V0 = [
      (V_plasma .* ( 1 - u.^2).^0.5 .* cos(phi))'
      (V_plasma .* ( 1 - u.^2).^0.5 .* sin(phi))'
      (V_plasma .* u)'
    ];
    
%    %delta function in plasma frame
    % V_plasma = (vecnorm(V_SW,2,1)) .* ones(1,sample_size); 
    % V0 = [
    %   V_plasma .* ( 1 - u'.^2).^0.5 .* cos(phi)'
    %   V_plasma .* ( 1 - u'.^2).^0.5 .* sin(phi)'
    %   V_plasma .* u'
    % ];
    
    
    % move into shock frame
    V0 = V0 + V_SW; %
    
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
