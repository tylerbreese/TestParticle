a = readtable('position_data.csv');
r = table2array(a);
b = readtable('velocity_data.csv');
v = table2array(b);
c = readtable('magfield_data.csv');
B = table2array(c);
e = readtable('distribution_data.csv');
En = table2array(e);
bins = 0:1:10;
En0 = En(:,4) ./ 1.6e-9;
Enf = En(:,5) ./ 1.6e-9;


R = 8.35e9;
figure()
scatter3(r(:,2)./R, r(:,3)./R, r(:,4)./R, [], r(:,1));
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
grid on;
% 
V = 400e5;
figure()
scatter3(v(:,2)./V, v(:,3)./V, v(:,4)./V, [], v(:,1));
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
grid on;

figure()
plot(r(:,2)./R,vecnorm(v(:,2:4),2,2)./V)

b0 = 0.05e-5;
figure()
hold on 
grid on
plot(r(:,2)./R,B(:,2)./b0)
plot(r(:,2)./R,B(:,3)./b0)
plot(r(:,2)./R,B(:,4)./b0)
plot(r(:,2)./R,B(:,5)./b0)

figure()
hold on
histogram(En0)
histogram(Enf)