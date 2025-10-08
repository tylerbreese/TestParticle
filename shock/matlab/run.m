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
global c e Z A q m
c = 3E10; %speed of light, m/s
e = 4.8E-10; %elementary charge, cgs
Z = 1; %atomic number
A = 1; %mass number
q = 1*e; %ion charge 
m = A * 1.67E-24; %g

U1 = vecnorm(V_SW,2,1);
B1 = vecnorm(B0,2,1);
Om = q*B1 ./ (m*c);
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
dt = 0.05 * (1/Om); %might need smaller timestep for shock frame
n = 1; t = 0.0;
Split = {};
for ii = 1:sample_size
    [B,U] = S.shock_field(B0, V_SW, x0, X(ii,:,1), An, kk, P, th, del);

    %while t < cycles*(1/Om)
    for n = 1:num_steps-1
        %integrate
        [X(ii,:,n+1),V(ii,:,n+1)] = integrate_boris(X(ii,:,n),V(ii,:,n),U,B,dt);
        %advance field
        [B,U] = S.shock_field(B0,V_SW,x0,X(ii,:,n+1),An,kk, P, th, del);
        Vmag(:,n+1) = sqrt( V(:,1,n+1).^2 + V(:,2,n+1).^2 + V(:,3,n+1).^2 );
        %cutoff = logspace(1,4,10);
        cutoff = [2,5,logspace(1,2,8)];
        [Split] = particle_split(X,V,Vmag,cutoff,Split);
        [Split] = particle_displace(Split,Rg);
        for jj = 1:length(Split)
            if isempty(Split{jj})
                [X_split,V_split] = integrate_boris(Split{jj}{1},Split{jj}{2},U,B,dt);
                Split{jj}{1} = X_split; Split{jj}{2} = V_split;
            else
                continue
            end
        end
        t = t + dt;
    end
end
toc
%%

function [X,V] = integrate_boris(Xn,Vn,U,B,dt)
    global q m c
    E(:,1) = (-1/c) .* (U(:,2).*B(:,3) - U(:,3).*B(:,2));
    E(:,2) = (-1/c) .* (U(:,3).*B(:,1) - U(:,1).*B(:,3));
    E(:,3) = (-1/c) .* (U(:,1).*B(:,2) - U(:,2).*B(:,1));
    % E = zeros(sample_size,3); %turn off E field
    %BORIS
    u = Vn;
    R = (q/m) .* E .* (dt./2);
    T = (q/(m*c)) .* B .* (dt./2);
    S = 2.*T ./ (1 + sum(T.*T,2));
    
    u = u + R; %first half step
    
    Vp(:,1) = u(:,1) + ( u(:,2).*T(:,3) - u(:,3).*T(:,2) );
    Vp(:,2) = u(:,2) + ( u(:,3).*T(:,1) - u(:,1).*T(:,3) ); %twirling
    Vp(:,3) = u(:,3) + ( u(:,1).*T(:,2) - u(:,2).*T(:,1) );
    
    Vq(:,1) = u(:,1) + ( Vp(:,2).*S(:,3) - Vp(:,3).*S(:,2) );
    Vq(:,2) = u(:,2) + ( Vp(:,3).*S(:,1) - Vp(:,1).*S(:,3) ); %twirling
    Vq(:,3) = u(:,3) + ( Vp(:,1).*S(:,2) - Vp(:,2).*S(:,1) );
    
    V = Vq + R; %next half step
    
    X(:,1) = Xn(:,1) + V(:,1) .* dt;
    X(:,2) = Xn(:,2) + V(:,2) .* dt;
    X(:,3) = Xn(:,3) + V(:,3) .* dt;
end

function [split] = particle_split(X,V,Vmag,cutoff,split)
    global m
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
