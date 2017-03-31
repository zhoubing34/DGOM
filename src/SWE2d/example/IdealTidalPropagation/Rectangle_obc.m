% read open boundary file
fp = fopen('Channel.edge', 'r');
Ns = fscanf(fp, '%d %d', 2); Ns = Ns(1);
data = fscanf(fp, '%d %d %d %d\n', [4, Ns]);
fclose(fp);

% obc vertex
surfID = data(4, :);
ind4 = find(surfID == 5); % east bc
vert1 = data(2, ind4);
vert2 = data(3, ind4);
vert = unique([vert1, vert2]);
Nv = numel(vert);

% obc data
T = 360;
w = 2*pi/T;
ftime = 10*T;
frequence = 12;
nt = 24*frequence;
t = linspace(0, ftime, nt);

fp = fopen('Rectangle.obc', 'w');
fprintf(fp, '%d %d\n', [Nv, 4]);
fprintf(fp, '%d ', vert);
fprintf(fp, '\n');
H = 40;
eta = 2;
phi = 0;
list = 1:Nv;
for i = 1:nt
    h = H + eta*cos(w*t(i) + phi);
    u = eta*sqrt(9.81/H)*cos(w*t(i) + phi);
    qx = u.*h;
    qy = 0;
    bot = 0;
    time = ones(1, Nv)*t(i);
    data = [list; time; ones(1, Nv)*h; ones(1, Nv)*qx; ones(1, Nv)*qy;...
        ones(1, Nv)*bot];
    fprintf(fp, '%d %20.16f %20.16f %20.16f %20.16f %20.16f\n', data);
end

fclose(fp);
