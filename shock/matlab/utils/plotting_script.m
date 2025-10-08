%%
% 
% h5create('output.h5','/position',size(X));
% h5write('output.h5','/position',X);
% h5create('output.h5','/velocity',size(V));
% h5write('output.h5','/velocity',V);
% h5create('output.h5','/gyro',[sample_size,2]);
% h5write('output.h5','/gyro',[r_g, w_g]);
% h5create('output.h5','/particle',[1,4]);
% h5write('output.h5','/particle',[m,q,Z,A]);
sample_size = 10;
X = h5read('output.h5','/position');
V = h5read('output.h5','/velocity');
gyroData = h5read('output.h5','/gyro');
particleData = h5read('output.h5','/particle');
m = particleData(1);
r_g = gyroData(:,1);
w_g = gyroData(:,2);


Vmag = vecnorm(V,2,2);

%% plotting
p = [0.4940 0.1840 0.5560]; %nice purple line
%pick = randi(sample_size); %grab any traj
pick = randsample(find(X(:,1,end)<0),1); %grab an ACR

tit1 = sprintf('V_0 = %2.1f km/s',Vmag(pick,1)/1e5);
r_gn = r_g(pick,1);
%r_gn = 1E12;
figure()
hold on
grid on
%plot(squeeze(X(pick,1,:))./1000,squeeze(X(pick,2,:))./1000,'.-','Color',p)
%plot3(squeeze(X(pick,1,:))./1000,squeeze(X(pick,2,:))./1000,squeeze(X(pick,3,:))./1000,'.-','Color',p)
plot3(squeeze(X(pick,1,1))./r_gn,squeeze(X(pick,2,1))./r_gn,squeeze(X(pick,3,1))./r_gn,'x-b') %plot start
plot3(squeeze(X(pick,1,end))./r_gn,squeeze(X(pick,2,end))./r_gn,squeeze(X(pick,3,end))./r_gn,'x-r') %plot end
plot3(squeeze(X(pick,1,:))./r_gn,squeeze(X(pick,2,:))./r_gn,squeeze(X(pick,3,:))./r_gn,'.-','Color',p)
title(tit1)
xlabel('X (r_{g})')
ylabel('Y (r_{g})')
zlabel('Z (r_{g})')
%%

U1 = vecnorm(V_SW,2,1);
B1 = vecnorm(B0,2,1);

particleData = h5read('output.h5','/particle');
m = particleData(1);
q = particleData(2);
Z = particleData(3);
A = particleData(4);

Om = q*B1 ./ (m*3e10);
Rg = U1 ./ Om;

E0 = (0.5*m).*Vmag(:,1).^2 ./ 1.6e-12;
Ef = (0.5*m).*Vmag(:,end).^2 ./ 1.6e-12;
posi = squeeze(X(:,1,1))./Rg;
posf = squeeze(X(:,1,end))./Rg;

figure()
hold on
grid on
scatter(E0./1000,posi,'o')
scatter(Ef./1000,posf,'x')
ylabel('X positions [r_{g}]')
xlabel('Energy [keV]');
title('Energy vs Final Position Shock 90deg to normal, shock at x=0, \tau=10 \Omega_{g}^{-1}');
legend('Initial Energy','Final Energy')

%%
%tau = linspace(1,t(pick),n)./t(pick);
%tau = linspace(1,n,n);
labx = sprintf('Time $(\\tau_{g})$');
laby = sprintf('Velocity (km/s)');
figure()
hold on
grid on
plot(Vmag(pick,:)./100000,'o-')
title(tit1)
xlabel(labx,'Interpreter','latex')
ylabel(laby)
ax = gca;
%ax.XScale = 'log';
ax.YScale = 'log';

% figure()
% hold on
% grid on
% histogram(reshape(Vmag(:,1),[],1),50)
% histogram(reshape(Vmag(:,end),[],1),50)
% title('Velocity Distribution')
% xlabel('cm/s')
% ax = gca;
% ax.XScale = 'log';
% %ax.YScale = 'log';
%%
[start,bins1] = histcounts(reshape(Vmag(:,1),[],1),50);
[finish,bins2] = histcounts(reshape(Vmag(:,end),[],1),50);
bins1 = (0.5*m).*(bins1).^2 ./ 1.6e-9;
bins2 = (0.5*m).*(bins2).^2 ./ 1.6e-9;
figure()
hold on
grid on
% histogram(reshape(Vmag(:,1),[],1),50)
% histogram(reshape(Vmag(:,end),[],1),50)
plot(bins1(2:end),start./sum(start))
plot(bins2(2:end),finish./sum(finish))
title('Velocity Distribution')
xlabel('keV')
ylabel('arbitrary units')
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';

%%
pick = randsample(find(X(:,1,end)<0),1); %grab an ACR

labx = sprintf('Position $(R/r_{g})$');
laby = sprintf('Velocity $(V/U_{SW})$');
vv = Vmag(pick,:)./vecnorm(V_SW,2,1); vv = vv(1:end-1)';
xx = squeeze(X(pick,1,:))./Rg; xx = xx(2:end-1);
% rr = [squeeze(X(pick,1,:)),squeeze(X(pick,2,:)),squeeze(X(pick,3,:))];
% rr = rr(:,1).^2 + rr(:,2).^2 + rr(:,3).^2;
% rr = (rr).^(1/2) ./ r_gn;

figure()
hold on
grid on
%plot(rr(1:end-1,1),vv,'.-')
plot(xx,vv,'.-')
title(tit1)
xlabel(labx,'Interpreter','latex')
ylabel(laby,'Interpreter','latex')
ax = gca;
%ax.XScale = 'log';
ax.YScale = 'log';

%% plotting
p = [0.4940 0.1840 0.5560]; %nice purple line
%pick = randi(sample_size); %grab any traj
pick = randsample(find(X(:,1,end)<0),1); %grab an ACR

tit1 = sprintf('V_0 = %2.1f km/s',Vmag(pick,1)/1e5);
r_gn = r_g(pick,1);
%r_gn = 1E12;
figure()
hold on
grid on
%plot(squeeze(X(pick,1,:))./1000,squeeze(X(pick,2,:))./1000,'.-','Color',p)
%plot3(squeeze(X(pick,1,:))./1000,squeeze(X(pick,2,:))./1000,squeeze(X(pick,3,:))./1000,'.-','Color',p)
plot3(squeeze(X(pick,1,1))./r_gn,squeeze(X(pick,2,1))./r_gn,squeeze(X(pick,3,1))./r_gn,'x-b') %plot start
plot3(squeeze(X(pick,1,end))./r_gn,squeeze(X(pick,2,end))./r_gn,squeeze(X(pick,3,end))./r_gn,'x-r') %plot end
plot3(squeeze(X(pick,1,:))./r_gn,squeeze(X(pick,2,:))./r_gn,squeeze(X(pick,3,:))./r_gn,'.-','Color',p)
title(tit1)
xlabel('X (r_{g})')
ylabel('Y (r_{g})')
zlabel('Z (r_{g})')

E0 = (0.5*m).*Vmag(:,1).^2 ./ 1.6e-12;
Ef = (0.5*m).*Vmag(:,end).^2 ./ 1.6e-12;
posi = squeeze(X(:,1,1))./r_g;
posf = squeeze(X(:,1,end))./r_g;

figure()
hold on
grid on
scatter(E0./1000,posi,'o')
scatter(Ef./1000,posf,'x')
ylabel('X positions [r_{g}]')
xlabel('Energy [keV]');
title('Energy vs Final Position Shock 80deg to normal, shock at x=0, \tau=100');
legend('Initial Energy','Final Energy')

%%
%tau = linspace(1,t(pick),n)./t(pick);
%tau = linspace(1,n,n);
labx = sprintf('Time $(\\tau_{g})$');
laby = sprintf('Velocity (km/s)');
figure()
hold on
grid on
plot(Vmag(pick,:)./100000,'o-')
title(tit1)
xlabel(labx,'Interpreter','latex')
ylabel(laby)
ax = gca;
%ax.XScale = 'log';
ax.YScale = 'log';

figure()
hold on
grid on
histogram(reshape(Vmag(:,1),[],1),50)
histogram(reshape(Vmag(:,end),[],1),50)
title('Velocity Distribution')
xlabel('cm/s')
ax = gca;
ax.XScale = 'log';
%ax.YScale = 'log';

%pick = randsample(find(X(:,1,end)<0),1); %grab an ACR
%pick = randsample(sample_size,1); %grab an ACR
labx = sprintf('Position $(x/r_{g})$');
laby = sprintf('Velocity $(v/U_{SW})$');
vv = Vmag(pick,:)./vecnorm(V_SW,2,1); vv = vv';
xx = squeeze(X(pick,1,:))./r_gn; xx = xx(2:end);
% rr = [squeeze(X(pick,1,:)),squeeze(X(pick,2,:)),squeeze(X(pick,3,:))];
% rr = rr(:,1).^2 + rr(:,2).^2 + rr(:,3).^2;
% rr = (rr).^(1/2) ./ r_gn;

figure()
hold on
grid on
%plot(rr(1:end-1,1),vv,'.-')
plot(xx,vv,'.-')
title(tit1)
xlabel(labx,'Interpreter','latex')
ylabel(laby,'Interpreter','latex')
ax = gca;
%ax.XScale = 'log';
ax.YScale = 'log';