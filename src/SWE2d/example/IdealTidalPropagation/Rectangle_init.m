fp = fopen('Channel.node', 'r');
Nv = fscanf(fp, '%d', 4); Nv = Nv(1);
data = fscanf(fp, '%d %f %f', [3, Nv]);
fclose(fp);

x = data(2, :);
y = data(3, :);

botLevel = zeros(size(x));

T = 360;
w = 2*pi/T;
H = 40;
g = 9.81;
c = sqrt(g*H);
k = w/c;
eta = .2;
phi = 0;
h = H + eta*cos(k*x);
u = eta*sqrt(g/H)*cos(k*x);
qx = u.*h;
qy = 0;
% write to initial file

data = zeros(5, numel(x));
data(1, :) = 1:Nv;
data(2, :) = h;
data(3, :) = qx;
data(4, :) = qy;
data(5, :) = botLevel;
fp = fopen('Rectangle.init', 'w');
fprintf(fp, '%d %d\n', Nv, 4);
fprintf(fp, '%d %20.16f %20.16f %20.16f %20.16f\n', data);
fclose(fp);
