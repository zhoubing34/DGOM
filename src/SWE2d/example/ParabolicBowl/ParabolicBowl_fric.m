fp = fopen('ParabolicBowl_Quad.node', 'r');
Nv = fscanf(fp, '%d', 4); Nv = Nv(1);
data = fscanf(fp, '%d %f %f', [3, Nv]);
fclose(fp);

x = data(2, :);
y = data(3, :);

Mann = 0;
data = zeros(2, Nv);
data(1, :) = 1:Nv;
data(2, :) = ones(1, Nv)*Mann;
fp = fopen('ParabolicBowl.fric', 'w');
fprintf(fp, '%d\n', Nv);
fprintf(fp, '%d %20.16f\n', data);
fclose(fp);