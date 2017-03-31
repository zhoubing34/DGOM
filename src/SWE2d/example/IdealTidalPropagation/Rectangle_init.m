fp = fopen('Channel.node', 'r');
Nv = fscanf(fp, '%d', 4); Nv = Nv(1);
data = fscanf(fp, '%d %f %f', [3, Nv]);
fclose(fp);

x = data(2, :);
y = data(3, :);

botLevel = zeros(size(x));
h = ones(size(x))*40;
% write to initial file

data = zeros(5, numel(x));
data(1, :) = 1:Nv;
data(2, :) = h;
data(5, :) = botLevel;
fp = fopen('Rectangle.init', 'w');
fprintf(fp, '%d %d\n', Nv, 4);
fprintf(fp, '%d %20.16f %20.16f %20.16f %20.16f\n', data);
fclose(fp);
