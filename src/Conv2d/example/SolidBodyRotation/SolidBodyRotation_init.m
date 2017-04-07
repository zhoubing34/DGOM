initialize_from_file = 'SolidBodyRotation.init';
node_file = 'SolidBodyRotation.node';

% read mesh
fp = fopen(node_file, 'r');
tmp = fscanf(fp, '%d %d %d %d', [4, 1]);
Nv = tmp(1);
data = fscanf(fp, '%d %f %f\n', [3, Nv]);
x = data(2,:);
y = data(3,:);
fclose(fp);

% obtain u,v and C
u = 0.5 - y;
v = x - 0.5;
r0 = 0.15;
c = zeros(size(x));
% slotted cylinder
x0 = 0.5;
y0 = 0.75;
r = sqrt((x-x0).^2 + (y-y0).^2)/r0;
ind = (r <= 1) & (( abs(x - x0)>0.025 ) | ( y>=0.85 ));
c(ind) = 1;

% cone
x0 = 0.5;
y0 = 0.25;
r = sqrt((x-x0).^2 + (y-y0).^2)/r0;
ind = (r <= 1);
c(ind) = 1 - r(ind);

% hump
x0 = 0.25;
y0 = 0.5;
r = sqrt((x-x0).^2 + (y-y0).^2)/r0;
ind = (r <= 1);
c(ind) = ( 1+cos(pi*r(ind)) )/4;

plot3(x, y, c, 'r.')

% write to file
data = [1:Nv; c; u; v];
fp = fopen(initialize_from_file, 'w');
fprintf(fp, '%d %d\n', Nv, 3);
fprintf(fp, '%d %20.16f %20.16f %20.16f\n', data);
fclose(fp);

