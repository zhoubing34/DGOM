function main
caseName = {'adv', 'diff'};
cellType = {'quad', 'tri'};
cellNum = [20, 40, 60, 80];
cellDeg = [1, 2, 3, 4];

caseid = [2, 1];
cellid = [3, 2];
dt_func = {@adv_dt, @diff_dt};

for icase = 1:numel(caseName)
    for itype = 1:numel(cellType)
        for inum = 1:numel(cellNum)
            for ideg = 1:numel(cellDeg)
                filename = [caseName{icase}, '_', cellType{itype}, '_',...
                    num2str(cellNum(inum)), '_', num2str(cellDeg(ideg))];
                fp = fopen([filename, '.inc'], 'w');
                fprintf(fp, '[============]DGOM: 2d convection problem\n');
                fprintf(fp, '[------------]case indicator (1 parameter)\n');
                fprintf(fp, '[------------]    1. case indicator  |-- 0. rotational convection\n');
                fprintf(fp, '[------------]                       |-- 1. advection-diffusion\n');
                fprintf(fp, '[------------]                       |-- 2. pure advection\n');
                fprintf(fp, '%d\n', caseid(icase));
                fprintf(fp, '[============]standard element info (2 parameters)\n');
                fprintf(fp, '[------------]    1. element type  |--- 2. triangle\n');
                fprintf(fp, '[------------]                     |--- 3. quadrilateral\n');
                fprintf(fp, '[------------]    2. order of polynomial;\n');
                fprintf(fp, '%d\n', cellid(itype));
                fprintf(fp, '%d\n', cellDeg(ideg));
                fprintf(fp, '[============]mesh info (3 parameters)\n');
                fprintf(fp, '[------------]    1. num of elements in x direction;\n');
                fprintf(fp, '[------------]    2. num of elements in y direction;\n');
                fprintf(fp, '%d\n', cellNum(inum));
                fprintf(fp, '%d\n', cellNum(inum));
                fprintf(fp, '[============]time info (3 parameters)\n');
                fprintf(fp, '[------------]    1. CFL number;\n');
                fprintf(fp, '[------------]    2. dt;\n');
                fprintf(fp, '[------------]    3. final time;\n');
                fprintf(fp, '%f\n', 0.1);
                dt = dt_func{icase}(cellNum(inum), cellDeg(ideg));
                fprintf(fp, '%f\n', dt);
                fprintf(fp, '%f\n', 2);
                fprintf(fp, '[============]parameters for advection-diffusion case\n');
                fprintf(fp, '[------------]    1. flow rate u for x direction;\n');
                fprintf(fp, '[------------]    2. flow rate v for x direction;\n');
                fprintf(fp, '[------------]    3. viscosity miu;\n');
                fprintf(fp, '%f\n', 0.5);
                fprintf(fp, '%f\n', 0.5);
                fprintf(fp, '%f\n', 0.01);
                fprintf(fp, '[============]output file\n');
                fprintf(fp, '[------------]    1. output file name;\n');
                fprintf(fp, '[------------]    2. output time interval;\n');
                fprintf(fp, '%s\n', filename);
                fprintf(fp, '0.001\n');
                fclose(fp);
            end% for
        end% for
    end% for
end% for

end

function dt = adv_dt(cellNum, deg)
spe = 0.5*sqrt(2);
len = 2/cellNum/(deg+1);
dt = len/spe;
end

function dt = diff_dt(cellNum, deg)
spe = 0.5*sqrt(2);
miu = 0.01;
len = 2/cellNum/(deg+1);
dt = min(len/spe, len^2/miu);
end


