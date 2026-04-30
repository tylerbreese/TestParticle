a = readtable('output/position_data.csv');
r = table2array(a);
b = readtable('output/velocity_data.csv');
v = table2array(b);
c = readtable('output/magfield_data.csv');
B = table2array(c);
e = readtable('output/distribution_data.csv');
En = table2array(e);
En0 = En(:,4);% ./ 1.6e-9;
Enf = En(:,5);% ./ 1.6e-9;

%%
T = r(:,1)*(0.05/6.78e6);
R = 1.0e5;
figure()
scatter3(r(:,2)./R, r(:,3)./R, r(:,4)./R, [], T);
xlabel('X-axis (km)');
ylabel('Y-axis (km)');
zlabel('Z-axis (km)');
grid on;
title('Position')
colorbar
c = colorbar;
c.Label.String = 'Time Elapsed (s)';
% 
V = 3e10;
figure()
scatter3(v(:,2)./V, v(:,3)./V, v(:,4)./V, [], T);
xlabel('X-axis (vx/c)');
ylabel('Y-axis (vy/c)');
zlabel('Z-axis (vz/c)');
title('Velocity')
grid on;
c = colorbar;
c.Label.String = 'Time Elapsed (s)';

figure()
scatter(r(:,2)./R,vecnorm(v(:,2:4),2,2)./V, [], T,'filled')
xlabel('Displacement (km)')
ylabel('Speed (v/c)')
c = colorbar;
c.Label.String = 'Time Elapsed (s)';
title('Phase Space')

b0 = 0.5;
figure()
hold on 
grid on
plot(r(:,2)./R,B(:,2)./b0)
plot(r(:,2)./R,B(:,3)./b0)
plot(r(:,2)./R,B(:,4)./b0)
plot(r(:,2)./R,B(:,5)./b0)
title('magnetic field')
%%
[n0,bins0] = histcounts(En0,25);
[n1,bins1] = histcounts(Enf,25);
n0(end+1) = 0.0; n1(end+1) = 0.0;
figure()
hold on
grid on
scatter(bins0,n0./sum(n0))
scatter(bins1,n1./sum(n1))
set(gca,'XScale','log','YScale','log');

%%