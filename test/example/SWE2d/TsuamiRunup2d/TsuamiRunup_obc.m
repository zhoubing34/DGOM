openind = 4;
obcfile = ['TsuamiRunup.obc',num2str(openind)];
edgefile = 'TsuamiRunup.edge';

% get vertex list
fp = fopen(edgefile, 'r');
tmp = fscanf(fp, '%d %d\n', [2, 1]);
nv = tmp(1);
SFToV = fscanf(fp, '%d %d %d %d', [4, nv]);
fclose(fp);

ind = (SFToV(4, :) == openind);
vert = SFToV([2,3], ind);
vert = unique(vert);
nv = numel(vert);

% get bottom elevation
botfile = 'TsuamiRunup.bot';
fp = fopen(botfile, 'r');
tmp = fscanf(fp, '%d\n', 1);
bot = fscanf(fp, '%d %lf', [2, tmp]);
bot = bot(2, :)';
fclose(fp);

% read open boundary condition
inputFile = 'Benchmark_2_input.txt';
fp = fopen(inputFile);
fgetl(fp);
inWave = fscanf(fp, '%e %e', [2, inf]);
fclose(fp);
time = inWave(1, :);
nt = numel(time);
eta = inWave(2, :);

% write obc file
fp = fopen(obcfile, 'w');
fprintf(fp, '%d %d\n', [nv, 3]);
fprintf(fp, '%d ', vert);
fprintf(fp, '\n');
for i = 1:nt
    tloc = time(i);
    etaloc = eta(i);
    hloc = etaloc - bot(vert);
    data = [1:nv; repmat(tloc, 1, nv); hloc'; zeros(1, nv); zeros(1, nv)];
    fprintf(fp, '%d %10.6e %20.16f %20.16f %20.16f\n', data);
end% for
fclose(fp);