fp = fopen('FlowOver3Bump.node', 'r');
Nv = fscanf(fp, '%d', 4); Nv = Nv(1);
data = fscanf(fp, '%d %f %f', [3, Nv]);
fclose(fp);

x = data(2, :);
y = data(3, :);

botLevel = zeros(size(x));

x0 = 30; y0 = 7.5; r0 = 7;
r2 = sqrt((x-x0).^2 + (y-y0).^2)./r0;
sk = r2 < 1.0;
botLevel(sk) = 1 - r2(sk);

x0 = 30; y0 = -7.5; r0 = 7;
r2 = sqrt((x-x0).^2 + (y-y0).^2)./r0;
sk = r2 < 1.0;
botLevel(sk) = 1 - r2(sk);

x0 = 47.5; y0 = 0; r0 = 8;
r2 = sqrt((x-x0).^2 + (y-y0).^2)./r0;
sk = r2 < 1;
botLevel(sk) = 2.8*(1 - r2(sk));
plot3(x, y, botLevel, '.')

% write to initial file
xb = 16;
h = zeros(size(x));
ind = x < xb;
h(ind) = 1.875;

data = zeros(5, numel(x));
data(1, :) = 1:Nv;
data(2, :) = h;
data(5, :) = botLevel;
fp = fopen('FlowOver3Bump.init', 'w');
fprintf(fp, '%d %d\n', Nv, 4);
fprintf(fp, '%d %20.16f %20.16f %20.16f %20.16f\n', data);
fclose(fp);
