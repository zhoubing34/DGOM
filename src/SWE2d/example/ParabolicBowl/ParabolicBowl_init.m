fp = fopen('ParabolicBowl_Quad.node', 'r');
Nv = fscanf(fp, '%d', 4); Nv = Nv(1);
data = fscanf(fp, '%d %f %f', [3, Nv]);
fclose(fp);

x = data(2, :);
y = data(3, :);

g = 9.81
alpha = 1.6*1e-7;
w     = sqrt(8*g*alpha);
T     = 2*pi/w;
X     = 1;
Y     = -0.41884;
r2    = (x.^2 + y.^2);
bot   = alpha*r2;
r2ext = (X+Y)/(alpha*(X^2 - Y^2));
h     = 1/(X+Y) + alpha*(Y^2 - X^2).*r2/(X+Y)^2;
h(r2>r2ext) = 0;
plot3(x, y, bot, 'b.'); hold on;
plot3(x, y, h+bot, 'r.');

% write to initial file

data = zeros(5, numel(x));
data(1, :) = 1:Nv;
data(2, :) = h;
data(5, :) = bot;
fp = fopen('ParabolicBowl.init', 'w');
fprintf(fp, '%d %d\n', Nv, 4);
fprintf(fp, '%d %20.16f %20.16f %20.16f %20.16f\n', data);
fclose(fp);
