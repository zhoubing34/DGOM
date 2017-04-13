function conv_rate
caseName = {'adv', 'diff'};
cellType = {'tri', 'quad'};
cellDeg = [1, 2, 3, 4];
cellNum = [20, 40, 60, 80];

num = numel(cellNum);

marker = {{'b--o', 'b--*', 'b--s', 'b--^'}, {'r-o', 'r-*', 'r-s', 'r-^'}};
markersize = 8;
markerface = {{'MarkerSize', markersize}, {'MarkerFaceColor', 'r', 'MarkerSize', markersize}};
fontsize = 16;
for icase = 1:numel(caseName)
    figure('color','w');
    for itype = 1:numel(cellType)
        for ideg = 1:numel(cellDeg)
            
            dx = zeros(num, 1);
            time = zeros(num, 1);
            L1 = zeros(num, 1);
            L2 = zeros(num, 1);
            Lf = zeros(num, 1);
            for inum = 1:num
                dx(inum) = 2/cellNum(inum);
                filename = [caseName{icase}, '_', cellType{itype}, '_', ...
                    num2str(cellNum(inum)), '_', num2str(cellDeg(ideg)), '.log'];
                [time(inum), L1(inum), L2(inum), Lf(inum)] = read_reasult(filename);
            end
            plot(dx, L2, marker{itype}{ideg}, markerface{itype}{:}); hold on;
        end
    end
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('$\Delta x$', 'Interpreter', 'Latex', 'FontSize', fontsize);
%     ylabel('$L_1$', 'Interpreter', 'Latex', 'FontSize', fontsize);
%     ylabel('$L_2$', 'Interpreter', 'Latex', 'FontSize', fontsize);
%     ylabel('$L_{\infty}$', 'Interpreter', 'Latex', 'FontSize', fontsize);
    ylabel('Ê±¼ä (s)', 'Interpreter', 'Latex', 'FontSize', fontsize);
end

end

function [time, L1, L2, Linf] = read_reasult(filename)
fp = fopen(filename, 'r');
% count line num
lineNum = 0;
while ~feof(fp)
    fgetl(fp);
    lineNum = lineNum + 1;
end
frewind(fp);
for i = 1:lineNum-2
    fgetl(fp);
end
% get elapsed time
str = fgetl(fp);
tid = find(str == ':');
time = str2double(str( (tid(2)+1):end ));
% get errors
str = fgetl(fp);
tid = find(str == ':');
eid = find(str == ',');
L1 = str2double(str( (tid(2)+1):(eid(2)-1) ));
L2 = str2double(str( (tid(3)+1):(eid(3)-1) ));
Linf = str2double(str( (tid(4)+1):end ));
fclose(fp);
end
