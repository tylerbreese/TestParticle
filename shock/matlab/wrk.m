pick = randi(sample_size);

x = squeeze(X(pick,1,:))./Rg;
y = squeeze(X(pick,2,:))./Rg;
z = squeeze(X(pick,3,:))./Rg;

vx   = squeeze(V(pick,1,:));
vy   = squeeze(V(pick,2,:));
vz   = squeeze(V(pick,3,:));
vmag = sqrt(vx.^2 + vy.^2 + vz.^2);

En = (0.5 * K.m) .* vmag.^2;
En = En ./ 1.6e-6; %in MeV

figure()
plot3(x,y,z)
figure()
plot(x,vmag);
figure()
histogram(En(:,end))
%%
tic
% --- constants ---
dt = 0.05 * (1/Om);
num_steps = cycles / dt;

% --- distribute particles among workers ---
spmd
    % Split X and V by particle index (dimension 1)
    Xd = codistributed(X);  % [num_particles x 3 x num_steps]
    Vd = codistributed(V);

    % Each worker gets its local chunk
    X_local = getLocalPart(Xd);
    V_local = getLocalPart(Vd);

    % Local parameters
    Split_local = {};
    [B,U] = S.shock_field_vec(B0, V_SW, x0, X_local(:,:,1), An, kk, P, th, del);

    % Local integration loop
    for n = 1:num_steps-1
        Xold = X_local(:,:,n);
        Vold = V_local(:,:,n);

        [Xnew,Vnew] = integrate_boris(Xold,Vold,U,B,dt,K);
        [B,U] = S.shock_field(B0,V_SW,x0,Xnew,An,kk,P,th,del);

        Vmag = sqrt( Vnew(:,1).^2 + Vnew(:,2).^2 + Vnew(:,3).^2 );

        X_local(:,:,n+1) = Xnew;
        V_local(:,:,n+1) = Vnew;

        cutoff = [2,5,logspace(1,2,8)];
        [Split_local] = particle_split(X_local,V_local,Vmag,cutoff,Split_local,K);
        [Split_local] = particle_displace(Split_local,Rg);

        for jj = 1:length(Split_local)
            if isempty(Split_local{jj})
                [X_split,V_split] = integrate_boris(Split_local{jj}{1},Split_local{jj}{2},U,B,dt,K);
                Split_local{jj}{1} = X_split;
                Split_local{jj}{2} = V_split;
            end
        end
    end

    % Combine results back into codistributed form
    Xd = codistributed.build(X_local, getCodistributor(Xd));
    Vd = codistributed.build(V_local, getCodistributor(Vd));
end

% --- Gather to client if needed ---
X = gather(Xd);
V = gather(Vd);
toc
%%
tic
numParticles = sample_size;
numWorkers = 4;  % or use numlabs inside spmd
chunk = ceil(numParticles / numWorkers);

spmd
    K.c = 3E10; %speed of light, m/s
    K.e = 4.8E-10; %elementary charge, cgs
    K.Z = 1; %atomic number
    K.A = 1; %mass number
    K.q = 1*K.e; %ion charge
    K.m = K.A * 1.67E-24; %g

    % Split up particle indices across workers
    labStart = (spmdIndex-1)*chunk + 1;
    labEnd = min(labStart + chunk - 1, numParticles);

    % Local particle slice
    X_local = X(labStart:labEnd,:,:);
    V_local = V(labStart:labEnd,:,:);

    % Local integration loop
    for n = 1:num_steps-1
        Xold = X_local(:,:,n);
        Vold = V_local(:,:,n);

        [Xnew,Vnew] = integrate_boris(Xold,Vold,U,B,dt,K);
        [B,U] = S.shock_field(B0,V_SW,x0,Xnew,An,kk,P,th,del);

        X_local(:,:,n+1) = Xnew;
        V_local(:,:,n+1) = Vnew;
    end
end

% Gather results
X = cat(1, X_local{:});
V = cat(1, V_local{:});
toc