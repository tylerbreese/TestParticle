tic
%% input
addpath('src');
% arg parser
if exist('input_file','var') && isfile(input_file)
    file = input_file;
else
    file = 'debug.in';
end
%matlab -batch "input_file='input.in'; run('script.m')"

fid = fopen(file, 'r');
%read in file
% fid         = fopen('debug.in','r');
formatSpec  = '%s %s %s\n';
inputdata   = textscan(fid,formatSpec);
run_name    = join([string(inputdata{1,3}{1}),'.h5'],'');
sample_size = str2double(inputdata{1,3}{2});
cycles      = str2double(inputdata{1,3}{3});
th          = str2double(inputdata{1,3}{4});
del         = str2double(inputdata{1,3}{5});
VSW         = str2double(inputdata{1,3}{6});
B_0         = str2double(inputdata{1,3}{7});
var         = str2double(inputdata{1,3}{8});

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
K.m = K.A * 9.11e-28; %g

U1 = vecnorm(V_SW,2,1);
B1 = vecnorm(B0,2,1);
Om = K.q*B1 ./ (K.m*K.c);
Rg = U1 ./ Om; 

%% initialize 
N = 200; %number of wave modes
S = ShockFields(B0,V_SW,[0, 0, 0],N);
P = InitParticles(V_SW, th, del, sample_size);

s = sqrt(var) .* B_0;
[An, kk] = S.init_turby(V_SW, Om, s, N);

%% integration
dt = 0.05 * (1/Om); %might need smaller timestep for shock frame
num_steps = ( cycles*(1/Om)/dt ) + 1; %number of steps

% X = zeros(sample_size,3,num_steps);
% V = zeros(sample_size,3,num_steps);
X = zeros(sample_size,3,3);
p = zeros(sample_size,3,3);
pmag = zeros(sample_size,3);

x0 = -20.*Rg .* ones(sample_size,1);
y0 = randi([-500,500],sample_size,1).*Rg;
z0 = randi([-500,500],sample_size,1).*Rg;

X(:,:,1) = [x0,y0,z0]; % start slighlt upstream of shock boundary
p(:, :, 1) = P.sampling(sample_size)';
pmag(:,1) = vecnorm(p(:,:,1),2,2);
pmag(:,3) = pmag(:,1);
X(:,:,3) = X(:,:,1);
p(:,:,3) = p(:,:,1);

%% integration
tic 
t = 0.0;
Split = {};
[B,U] = S.shock_field_vec(B0, V_SW, 0.0, X(:,:,1), An, kk, P, th, del);

%while t < cycles*(1/Om)
for n = 1:num_steps-1
    %integrate
    % Xold = X(:,:,n);
    % Vold = V(:,:,n);
    Xold = X(:,:,1);
    Vold = p(:,:,1);
    
    [Xnew,Vnew] = integrate_boris(Xold,Vold,U,B,dt,K);
    %advance field
    [B,U] = S.shock_field_multipole(B0,V_SW,0.0,Xnew,An,kk, P, th, del);
    pmag(:,2) = sqrt( Vnew(:,1).^2 + Vnew(:,2).^2 + Vnew(:,3).^2 );
    % X(:,:,n+1) = Xnew;
    % V(:,:,n+1) = Vnew;
    X(:,:,2) = Xnew;
    p(:,:,2) = Vnew;
    X(:,:,1) = X(:,:,2);
    p(:,:,1) = p(:,:,2);
    cutoff = [2,5,logspace(1,2,8)];
    % [Split] = particle_split(X,p,pmag,cutoff,Split,K);
    % [Split] = particle_displace(Split,Rg);
    % for jj = 1:length(Split)
    %     if isempty(Split{jj})
    %         [X_split,V_split] = integrate_boris(Split{jj}{1},Split{jj}{2},U,B,dt,K);
    %         Split{jj}{1} = X_split; Split{jj}{2} = V_split;
    %     else
    %         continue
    %     end
    % end
    t = t + dt;
end
toc
%% sim end save data
% qqx = [];
% for ii = 1:length(Split)
%     qqx = cat(1,qqx,Split{ii}{1});
% end
% qqp = [];
% for ii = 1:length(Split)
%     qqp = cat(1,qqp,Split{ii}{2});
% end
% qqf = [];
% for ii = 1:length(Split)
%     qqf = cat(1,qqf,Split{ii}{3});
% end
out_folder = "multi";
h5out = join([out_folder, '/', run_name],'');
h5create(h5out,'/position',size(X));
h5write(h5out,'/position',X);
h5create(h5out,'/velocity',size(p));
h5write(h5out,'/velocity',p);
h5create(h5out,'/gyro',[1,2]);
h5write(h5out,'/gyro',[Rg, Om]);
h5create(h5out,'/particle',[1,4]);
h5write(h5out,'/particle',[K.m,K.q,K.Z,K.A]);
% h5create(h5out,'/split_x',size(qqx));
% h5write(h5out,'/split_x',qqx);
% h5create(h5out,'/split_p',size(qqp));
% h5write(h5out,'/split_p',qqp);
% h5create(h5out,'/split_flag',size(qqf));
% h5write(h5out,'/split_flag',qqf);
%%


function [X,p] = integrate_boris(Xn,pn,U,B,dt,K)
    q = K.q;
    m = K.m;
    c = K.c;
    gamma = (1 + ((vecnorm(pn,2,2)./(m*c))).^2).^(0.5);
   
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

function [split] = particle_split(X,p,pmag,cutoff,split,K)
    m = K.m;
    c = 3e10;
    En = (0.5*m).*(pmag(:,2).^2);
    E0 = (0.5*m).*(pmag(:,3).^2);
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
            split{ii}{1} = X_split; split{ii}{2} = V_split;
        end
    end
end
