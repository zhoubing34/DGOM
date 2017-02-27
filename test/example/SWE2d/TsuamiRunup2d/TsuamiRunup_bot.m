% read input topography
bathyfile = 'Benchmark_2_Bathymetry.txt';
fp = fopen(bathyfile);
fgetl(fp);
data = fscanf(fp, '%e %e %e\n', [3, inf]);
fclose(fp);
x_e = data(1, :)';
y_e = data(2, :)';
h_e = data(3, :)';

% read mesh
nodefile = 'TsuamiRunup.node';
fp = fopen(nodefile);
tmp = fscanf(fp, '%d %d %d %d', [4, 1]);
nv = tmp(1);
data = fscanf(fp, '%d %lf %lf', [3, inf]);
fclose(fp);
vx = data(2, :)';
vy = data(3, :)';

% create bot file
% interpolate
interp = TriScatteredInterp(x_e,y_e,-h_e,'linear');
bot = interp(vx, vy);
botfile = 'TsuamiRunup.bot';
fp = fopen(botfile, 'w');
fprintf(fp, '%d\n', nv);
fprintf(fp, '%d %20.16e\n', [1:nv; bot']);
fclose(fp);
plot3(vx, vy, bot, 'r.');