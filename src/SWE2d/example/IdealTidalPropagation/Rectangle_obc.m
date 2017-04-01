% read open boundary file
fp = fopen('Channel.edge', 'r');
Ns = fscanf(fp, '%d %d', 2); Ns = Ns(1);
data = fscanf(fp, '%d %d %d %d\n', [4, Ns]);
fclose(fp);

% obc vertex
surfID = data(4, :);
ind4 = find(surfID >= 4); % east bc
vert1 = data(2, ind4);
vert2 = data(3, ind4);
vert = unique([vert1, vert2]);
Nv = numel(vert);

% coordinate
fp = fopen('Channel.node', 'r');
Nvert = fscanf(fp, '%d', 4); Nvert = Nvert(1);
data = fscanf(fp, '%d %f %f', [3, Nvert]);
fclose(fp);
x = data(2, :);

% obc data
T = 360;
w = 2*pi/T;
ftime = 20*T;
frequence = 24;
nt = 10*frequence;
t = linspace(0, ftime, nt);

fp = fopen('Rectangle.obc', 'w');
fprintf(fp, '%d %d\n', [Nv, 4]);
fprintf(fp, '%d ', vert);
fprintf(fp, '\n');
H = 40;
g = 9.81;
c = sqrt(g*H);
k = w/c;
eta = .2;
phi = 0;
list = 1:Nv;
for i = 1:nt
    h = H + eta*cos(k*x(vert) - w*t(i));
    u = eta*sqrt(g/H)*cos(k*x(vert) - w*t(i));
    qx = u.*h;
    qy = 0;
    bot = 0;
    time = ones(1, Nv)*t(i);
    data = [list; time; h; qx; ones(1, Nv)*qy;...
        ones(1, Nv)*bot];
    fprintf(fp, '%d %20.16f %20.16f %20.16f %20.16f %20.16f\n', data);
end

fclose(fp);
