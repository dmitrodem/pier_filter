close all;
clear;
clc;

do_run = 0;
physical_constants;

unit = 1e-3;


f_0 = 35.0e9;
fbw = 0.1;
f_start = f_0*(1-2*fbw);
f_end   = f_0*(1+2*fbw);
lambda0 = c0/f_0/unit;

TE_mode = 'TE10';

FDTD = InitFDTD('NrTS', 1e7, 'OverSampling', 5, 'EndCriteria', 1e-10);
FDTD = SetupMPI(FDTD,'SplitN_X',2);
BC = {'PML_8', 'PML_8', 'PEC', 'PEC', 'PEC', 'PEC'};
FDTD = SetBoundaryCond(FDTD, BC);
FDTD = SetGaussExcite(FDTD, 0.5*(f_start + f_end), 0.5*(f_end - f_start));

mesh_res = lambda0 ./ [500, 500, 500];

gaps = [4.56, 3.59, 3.39, 3.36, 3.39, 3.59, 4.56];
lengths = [4.68, 5.50, 5.63, 5.63, 5.50, 4.68];
wg_width = 7.112;
wg_height = wg_width/2;
wg_length = 5.0;
thickness = 2.0;

CSX = InitCSX();

refine_y = horzcat(gaps/2, -gaps/2, wg_width/2, -wg_width/2, 0);
refine_x = [0]

CSX = AddMetal(CSX, 'PEC');

xpos = wg_length;
for i = 1:size(gaps, 2)
  start = [xpos,           gaps(i)/2,  0];
  stop  = [xpos+thickness, wg_width/2, wg_height];
  CSX = AddBox(CSX, 'PEC', 1, start, stop);
  start = [xpos,           -gaps(i)/2,  0];
  stop  = [xpos+thickness, -wg_width/2, wg_height];
  CSX = AddBox(CSX, 'PEC', 1, start, stop);
  refine_x(end + 1) = xpos;
  refine_x(end + 1) = xpos + thickness;
  if (i < size(gaps, 2))
    xpos  = xpos + lengths(i);
  end
end
refine_x(end + 1) = refine_x(end) + wg_length;

mesh.x = SmoothMeshLines(refine_x,       mesh_res(1));
mesh.y = SmoothMeshLines(refine_y,       mesh_res(2));
mesh.z = SmoothMeshLines([0, wg_height], mesh_res(3));

CSX = DefineRectGrid(CSX, unit, mesh);

start = [mesh.x(10), mesh.y(1),   mesh.z(1)];
stop  = [mesh.x(14), mesh.y(end), mesh.z(end)];
[CSX, port{1}] = AddRectWaveGuidePort(CSX, 0, 1, start, stop, 'x', wg_width*unit, wg_height*unit, TE_mode, 1);

start = [mesh.x(end-13), mesh.y(1),   mesh.z(1)];
stop  = [mesh.x(end-14), mesh.y(end), mesh.z(end)];
[CSX, port{2}] = AddRectWaveGuidePort(CSX, 0, 2, start, stop, 'x', wg_width*unit, wg_height*unit, TE_mode);

%CSX = AddDump(CSX,'Et','FileType',1,'SubSampling','2,2,2');
%start = [mesh.x(1)   mesh.y(1)   mesh.z(1)];
%stop  = [mesh.x(end) mesh.y(end) mesh.z(end)];
%CSX = AddBox(CSX, 'Et', 0, start, stop);

SimPath = 'wg_run';
SimCSX = 'wg.xml';

if (do_run == 1)
  [status, message, messageid] = rmdir(SimPath, 's');
  [status, message, messageid] = mkdir(SimPath);

  WriteOpenEMS([SimPath, '/', SimCSX], FDTD, CSX);
  %CSXGeomPlot([SimPath, '/', SimCSX]);
  %RunOpenEMS(SimPath, SimCSX);
end

