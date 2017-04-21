function SectionProfile(ratio)
% Compare the results with exact solution on section x=0;
g     = 9.81;
alpha = 1.6*1e-7;
w     = sqrt(8*g*alpha);
T     = 2*pi/w;
delta = 0.5;
rmin  = -4000+delta; 
rmax  =  4000-delta;

%% Get points coordinate
ne    = 100;   % number of exact solution
np    = 50;    % number of interpolated solutions
xe    = zeros(ne, 1);
ye    = linspace(rmin, rmax, ne)';
re    = (xe.^2 + ye.^2);
be    = alpha*re;

xp    = zeros(np, 1);  
yp    = linspace(rmin, rmax, np)';
rp    = (xp.^2 + yp.^2);
bp    = alpha*rp;

%% Set output file
file0 = 'ParabolicBowl_Quad.0-2.nc';
file1 = 'ParabolicBowl_Quad.1-2.nc';

%% spicific time
time0 = ncread(file0, 'time');
nt = numel(time0); % number of time step
tStep = round((nt-1)*ratio + 1);
% timeFrac = (1/4:1/4:3/4);
time     = time0(tStep);
timeStr  = cell(numel(tStep), 1);
for i = 1:numel(tStep)
    timeStr{i} = ['t = ',num2str( time(i)/T ),'T',];
end

%% draw figure
lineWidth = 2;
markerSize = 5;
fontSize = 12;

% get nodes coordinate
x0 = double(ncread(file0, 'x')); y0 = double(ncread(file0, 'y'));
x1 = double(ncread(file1, 'x')); y1 = double(ncread(file1, 'y'));

for ist = 1:numel(time)
    % get result
    h0 = ncread(file0, 'h', [1,1,tStep(ist)], [inf,inf,1]);
    h1 = ncread(file1, 'h', [1,1,tStep(ist)], [inf,inf,1]);
    q0 = ncread(file0, 'qy', [1,1,tStep(ist)], [inf,inf,1]);
    q1 = ncread(file1, 'qy', [1,1,tStep(ist)], [inf,inf,1]);

    % get exact solution
    [h_ext, ~, qy_ext] = ExtSolution(xe, ye, time(ist) );

    % interpolation
    intp0 = TriScatteredInterp(x0(:), y0(:), double(h0(:)));
    intp1 = TriScatteredInterp(x1(:), y1(:), double(h1(:)));
    h_inp0 = intp0(xp(:), yp(:));
    h_inp1 = intp1(xp(:), yp(:));
    [h_inp0, h_inp1] = SetCommonData(h_inp0, h_inp1);
    intp0 = TriScatteredInterp(x0(:), y0(:), double(q0(:)));
    intp1 = TriScatteredInterp(x1(:), y1(:), double(q1(:)));
    q_inp0 = intp0(xp(:), yp(:));
    q_inp1 = intp1(xp(:), yp(:));
    [q_inp0, q_inp1] = SetCommonData(q_inp0, q_inp1);
    
    % plot
    figure('Color', 'w'); subplot(2,1,1);
    plot(ye, h_ext+be, 'k-', 'LineWidth',lineWidth); hold on;
    plot(yp, h_inp0+bp, 'ro', 'MarkerSize', markerSize, ...
        'LineWidth',lineWidth, 'MarkerFaceColor', 'r');
    plot(yp, h_inp1+bp, 'bo', 'MarkerSize', markerSize, ...
        'LineWidth',lineWidth, 'MarkerFaceColor', 'b');

    ylabel('$\eta (m)$', 'Interpreter', 'latex', 'FontSize', fontSize);
    xlabel('$y (m)$', 'Interpreter', 'latex', 'FontSize', fontSize);
    legend({'Exact', 'proc:0', 'proc:1'}, 'box', 'off','FontSize', fontSize);
    plot(ye, be, 'k', 'LineWidth',lineWidth);

    subplot(2,1,2);
    plot(ye, qy_ext, 'k-', 'LineWidth',lineWidth); hold on;
    plot(yp, q_inp0, 'ro', 'MarkerSize', markerSize, ...
        'LineWidth',lineWidth, 'MarkerFaceColor', 'r');
    plot(yp, q_inp1, 'bo', 'MarkerSize', markerSize, ...
        'LineWidth',lineWidth, 'MarkerFaceColor', 'b');
    ylim([-1.25, 1.25]);
    ylabel('$q_y (m^2/s)$', 'Interpreter', 'latex', 'FontSize', fontSize);
    xlabel('$y (m)$', 'Interpreter', 'latex', 'FontSize', fontSize);
end% for
end% func

function [h_inp0, h_inp1] = SetCommonData(h_inp0, h_inp1)
common_flag = ~(isnan(h_inp0) | isnan(h_inp1));
h_common = max(h_inp0(common_flag), h_inp1(common_flag));
h_inp0(common_flag) = h_common;
h_inp1(common_flag) = h_common;
end