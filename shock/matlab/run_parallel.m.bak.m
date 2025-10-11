tic
%% input
%read in file
fid         = fopen('debug.in','r');
formatSpec  = '%s %s %s\n';
inputdata   = textscan(fid,formatSpec,7);
run_name    = join([string(inputdata{1,3}{1}),'_s0p5','.h5'],'');
sample_size = str2double(inputdata{1,3}{2});
cycles      = str2double(inputdata{1,3}{3});
th          = str2double(inputdata{1,3}{4});
del         = str2double(inputdata{1,3}{5});
VSW         = str2double(inputdata{1,3}{6});
B_0         = str2double(inputdata{1,3}{7});

% get initial cond
V_SW        = VSW .* [cos(del); 0.0; sin(del)]; %solar wind speed, cm/s
B0          = B_0 .* [cos(th); 0.0; sin(th)]; %magnetic field in Tesla B(r) = B0/R[AU]

%% constants
% global c e Z A q m
K.c = 3E10; %speed of light, m/s
K.e = 4.8E-10; %elementary charge, cgs
K.Z = 1; %atomic number
K.A = 1; %mass number
K.q = 1*K.e; %ion charge 
K.m = K.A * 1.67E-24; %g

U1 = vecnorm(V_SW,2,1);
B1 = vecnorm(B0,2,1);
Om = K.q*B1 ./ (K.m*K.c);
Rg = U1 ./ Om;
num_steps = cycles*(1/Om)*20+1; %number of steps 

%% initialize 
N = 200; %number of wave modes
S = ShockFields(B0,V_SW,[0, 0, 0],N);
P = InitParticles(V_SW, th, del, sample_size);

s = sqrt(0.5) .* B_0;
[An, kk] = S.init_turby(V_SW, Om, s, N);

X = zeros(sample_size,3,num_steps);
V = zeros(sample_size,3,num_steps);
Vmag = zeros(sample_size,num_steps);

x0 = -20.*Rg .* ones(sample_size,1);
y0 = randi([-500,500],sample_size,1).*Rg;
z0 = randi([-500,500],sample_size,1).*Rg;

X(:,:,1) = [x0,y0,z0]; % start slighlt upstream of shock boundary
V(:, :, 1) = P.sampling(sample_size, V_SW);
Vmag(:,1) = vecnorm(V(:,:,1),2,2);
%% integration
dt = 0.05 * (1/Om); % timestep

% Preallocate arrays
X_out = zeros(sample_size,3,num_steps);
V_out = zeros(sample_size,3,num_steps);

parfor (ii = 1:sample_size,8)
    
    % Initialize local variables
    t_local = 0;
    Split_local = {};
    X_local = zeros(3,num_steps);
    V_local = zeros(3,num_steps);

    % Initial conditions
    X_local(:,1) = X(ii,:,1).';
    V_local(:,1) = V(ii,:,1).';
    
    % Compute initial field
    [B,U] = S.shock_field(B0, V_SW, x0, X_local(:,1), An, kk, P, th, del);

    for n = 1:num_steps-1
        % Integrate particle positions/velocities
        [X_local(:,n+1), V_local(:,n+1)] = integrate_boris(X_local(:,n), V_local(:,n), U, B, dt, K);
        
        % Update fields at current position
        [B,U] = S.shock_field(B0, V_SW, x0, X_local(:,n+1), An, kk, P, th, del);
        
        % Compute speed magnitude
        Vmag = sqrt(sum(V_local(:,n+1).^2));
        
        % % Particle splitting/displacement
        % cutoff = [2,5,logspace(1,2,8)];
        % Split_local = particle_split(X_local, V_local, Vmag, cutoff, Split_local);
        % Split_local = particle_displace(Split_local, Rg);
        % 
        % % Integrate split particles
        % for jj = 1:length(Split_local)
        %     if ~isempty(Split_local{jj})
        %         [X_split, V_split] = integrate_boris(Split_local{jj}{1}, Split_local{jj}{2}, U, B, dt);
        %         Split_local{jj}{1} = X_split;
        %         Split_local{jj}{2} = V_split;
        %     end
        % end

    end

    % Save local results
    X_out(ii,:,:) = X_local;
    V_out(ii,:,:) = V_local;
    % Split_out{ii} = Split_local;
end
X = X_out;
V = V_out;
% Split = Split_out;

toc
%%

function [X,V] = integrate_boris(Xn,Vn,U,B,dt,K)
    q = K.q;
    m = K.m;
    c = K.c;
    % E(1) = (-1./c) .* (U(2).*B(3) - U(3).*B(2));
    % E(2) = (-1./c) .* (U(3).*B(1) - U(1).*B(3));
    % E(3) = (-1./c) .* (U(1).*B(2) - U(2).*B(1));
    E(1) = (c).^(-1) .* (U(2).*B(3) - U(3).*B(2));
    E(2) = (c).^(-1) .* (U(3).*B(1) - U(1).*B(3));
    E(3) = (c).^(-1) .* (U(1).*B(2) - U(2).*B(1));
    % E = zeros(sample_size,3); %turn off E field
    %BORIS
    u = Vn;

    R(1) = (q/m) .* E(1) .* (dt./2);
    R(2) = (q/m) .* E(2) .* (dt./2);
    R(3) = (q/m) .* E(3) .* (dt./2);

    T(1) = (q/(m*c)) .* B(1) .* (dt./2);
    T(2) = (q/(m*c)) .* B(2) .* (dt./2);
    T(3) = (q/(m*c)) .* B(3) .* (dt./2);

    S(1) = 2.*T(1) ./ (1 + sum(T.*T,2));
    S(2) = 2.*T(2) ./ (1 + sum(T.*T,2));
    S(3) = 2.*T(3) ./ (1 + sum(T.*T,2));
    
    u(1) = u(1) + R(1); %first half step
    u(2) = u(2) + R(2); %first half step
    u(3) = u(3) + R(3); %first half step
    
    Vp(1) = u(1) + ( u(2).*T(3) - u(3).*T(2) );
    Vp(2) = u(2) + ( u(3).*T(1) - u(1).*T(3) ); %twirling
    Vp(3) = u(3) + ( u(1).*T(2) - u(2).*T(1) );
    
    Vq(1) = u(1) + ( Vp(2).*S(3) - Vp(3).*S(2) );
    Vq(2) = u(2) + ( Vp(3).*S(1) - Vp(1).*S(3) ); %twirling
    Vq(3) = u(3) + ( Vp(1).*S(2) - Vp(2).*S(1) );

    V(1) = Vq(1) + R(1); %next half step
    V(2) = Vq(2) + R(2); %next half step
    V(3) = Vq(3) + R(3); %next half step
    
    X(1) = Xn(1) + V(1) .* dt;
    X(2) = Xn(2) + V(2) .* dt;
    X(3) = Xn(3) + V(3) .* dt;
end

function [split] = particle_split(X,V,Vmag,cutoff,split,K)
    m = K.m;
    En = (0.5*m).*(Vmag(:,end).^2);
    E0 = (0.5*m).*(Vmag(:,1).^2);
    ratio = En ./ E0;
    flag = ratio>cutoff;
    if isempty(split)
        for ii = 1:length(cutoff)
            split{ii} = {X(flag(:,ii),:,2), V(flag(:,ii),:,2), flag(:,ii)};
        end
    else
        temp_flag = flag;
        for ii = 1:length(split)
            flag(:,ii) = split{ii}{3};
        end
        flag = flag | temp_flag;
        flag(flag>1) = 1;
        for ii = 1:length(cutoff)
            ssplit{ii} = {X(flag(:,ii),:,2), V(flag(:,ii),:,2), flag(:,ii)};
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
            split{ii}{1} = X_split; split{ii}{2} = V_split;
        end
    end
end
