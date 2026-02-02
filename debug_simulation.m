% Debug Script for solvStitcher (Simple)
clear; clc;

% 1. Setup Same as run_simulation.m
D = 4.7625e-3;
dm = 20e-3;
Z = 14;

trulla = BallBearing();

fi = 0.52; fo = 0.52;
ri = fi*D; ro = fo*D;
alpha_free = 15;

% Corrected call
trulla = setGeometry(trulla, D, dm, alpha_free, ri, ro, Z, 0, 0);

rho_ball = 3200;
E_ball = 3.10e11; nu_ball = 0.26;
E_ring = 2.08e11; nu_ring = 0.3;

trulla = setPhysical(trulla, rho_ball, E_ball, E_ring, nu_ball, nu_ring);

trulla.lubrication.eta_0 = 0.05;
trulla.lubrication.alpha = 2.2e-8;
trulla.lubrication.rho = 900;

load_case = [100, 0, 0, 0, 0, 20000];

% 2. Internal Setup from calcDisplacement.m
bea_geo = trulla.geometry;
bea_phys = trulla.physical;
bea_load = makeLoad(load_case);

% 3. Initialize
init = solvInit(bea_geo.z);

% 4. Call solvStitcher Directy
disp('Calling solvStitcher directly...');
[res, jac] = solvStitcher(trulla, bea_geo, bea_phys, bea_load, 'outer', init);

disp('solvStitcher completed without error.');
disp('Result Norm:');
disp(norm(res));
