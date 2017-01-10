%% flow field
[x,y] = meshgrid(-1:0.1:1,-1:0.1:1);
w = 5*pi/6;
u = -w.*y;
v = w.*x;
z = x.*0;

figure; 
mesh(x,y,z); view(0, 90)
hold on;
q = quiver(x,y,u,v);
quiverkey(q, 0, 1.2, 5, 'm/s', 'Color', 'b', 'LabelDistance', 0.03)
grid on;
axis equal;
box on;
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);

%% initial value
t = linspace(-1, 1, 100);
[x,y] = meshgrid(t,t);
sigma = 125*1e3/33^2;
x0 = 0; y0 = 3/5;
z = exp(-sigma.*((x-x0).^2+(y-y0).^2));
figure; 
mesh(x,y,z);
hold on;
grid on;
axis equal;
box on;