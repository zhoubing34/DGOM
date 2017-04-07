folder = 'obcE';
file0 = [folder,'/swe2d_channel.0-2.nc'];
file1 = [folder,'/swe2d_channel.1-2.nc'];
x0 = ncread(file0, 'x'); y0 = ncread(file0, 'y');
x1 = ncread(file1, 'x'); y1 = ncread(file1, 'y');
x = x1; y = y1; file = file1;
[n, k] = find(  (abs(x - 2e4)<1e-4 & abs(y - 0)<1e-4) );
k0 = k(1);
n0 = n(1);
xe = x1(n0, k0);
ye = y1(n0, k0);
time = ncread(file0, 'time');
nt = numel(time);
h0 = zeros(nt, 1);
q0 = zeros(nt, 1);
for t = 1:nt
    h = ncread(file, 'h', [1,1,t], [inf,inf,1]);
    qx = ncread(file, 'qx', [1,1,t], [inf,inf,1]);
    h0(t) = h(n0, k0);
    q0(t) = qx(n0, k0);
end
% exact value
H = 40;
g = 9.81;
c = sqrt(g*H);
T = 360;
w = 2*pi/T;
k = w/c;
eta = .2;
h_ext = H + eta*cos(k*xe - w*time);
u_ext = eta*sqrt(g/H)*cos(k*xe - w*time);
q_ext = h_ext.*u_ext;
figure('color', 'w');
subplot(2,2,1); xlabel('time (s)'), ylabel('h'); grid on; hold on;
plot(time, h_ext, 'k.');
plot(time, h0, 'r.-');
subplot(2,2,3); xlabel('time (s)'), ylabel('qx'); grid on; hold on;
plot(time, q_ext, 'k.');
plot(time, q0, 'b.-');
%% closed bay
file0 = [folder,'/swe2d_closedbay.0-2.nc'];
file1 = [folder,'/swe2d_closedbay.1-2.nc'];
file = file1;
time = ncread(file0, 'time');
nt = numel(time);
h0 = zeros(nt, 1);
q0 = zeros(nt, 1);
for t = 1:nt
    h = ncread(file, 'h', [1,1,t], [inf,inf,1]);
    qx = ncread(file, 'qx', [1,1,t], [inf,inf,1]);
    h0(t) = h(n0, k0);
    q0(t) = qx(n0, k0);
end
h_ext = H + eta*(cos(k*xe - w*time) + cos(-k*xe + w*time));
% u_ext = eta*sqrt(g/H)*(cos(k*xe - w*time) + cos(-k*xe + w*time));
q_ext = zeros(size(time));
subplot(2,2,2); xlabel('time (s)'), ylabel('h'); grid on; hold on;
plot(time, h_ext, 'k.');
plot(time, h0, 'r.-');
subplot(2,2,4); xlabel('time (s)'), ylabel('qx'); grid on; hold on;
plot(time, q_ext, 'k.');
plot(time, q0, 'b.-');