function [resStruct, solVec, dispArray]=calcDisplacement(obj,...
    loadcase, rc_assumption, init_conditions)
%Main method to run a simulation and obtain a displacement.
%Retruns a formated struct (resStruct) of the results as well
%as raw data. Extracting the required informaiton from
%resStruct is the most user friendly way of using this
%function.
%solVec can directly be used as init_conditions for a similar
%load case

%loadcase contains the particular load applied for the
%simulation and must be an array of the following shape:
%[Axial load, radial Y force, radial Z force, angular error
%around the y-axis, angular error around the z-axis, rotational
%speed in RPM]

rc_scenario=rc_assumption;
if (obj.geometry.set==0)
    error('Bearing geometry not set. Use method setGeometry.');
end
if (obj.physical.set==0)
    error(['Bearing physical properties not set. Use method'...
        ' setPhysical.']);
end

bea_geo=obj.geometry;
bea_phys=obj.physical;

bea_load=makeLoad(loadcase);

if isempty(init_conditions)
    init=solvInit(bea_geo.z);
else
    init=init_conditions;
end

options = optimoptions('fsolve','Display','iter', ...
    'MaxFunEvals', 1e9,...
    'MaxIter', 1e9, ...
    'Jacobian', 'off',...
    'FunValCheck','on',...
    'TolFun',1e-18,...
    'TolX', 1e-18, ...
    'Algorithm','trust-region-dogleg');
% 'JacobPattern',jPat(bea_geo.z),... % only used by trust-region reflective
% levenberg-marquardt
% trust-region-reflective
% trust-region-dogleg
%'FinDiffRelStep', 1e-9,...
%  'FinDiffType', 'central',...

loadSolv=@(val)solvStitcher(obj, bea_geo, bea_phys, bea_load,...
    rc_scenario, val);

[solVec,fval,exitflag,output] = ...
    fsolve(loadSolv,init,options);

resStruct = genSolvStruct(solVec, bea_geo);

[resStruct,rc_inner] = addSigmaMaxAlpha(obj, resStruct, bea_geo,...
    bea_phys, bea_load);

%checks if the race control assumption was correct. If it was
%not, the solver is re-run. This only checks for one ball
%(with median stress) and decides the race control scenario for
%the entire bearing. It is a know limitation of this class,
%that race control is not determined for every ball separately.

if rc_inner
    if strcmp(rc_scenario,'outer')
        display('Inner race control detected. Re-running solver.');
        rc_scenario='inner';
        loadSolv=@(val)solvStitcher(obj, bea_geo, bea_phys, bea_load,...
            rc_scenario, val);
        [solVec,fval,exitflag,output]=fsolve(loadSolv,init,...
            options);
        resStruct = addSigmaMaxAlpha(obj, resStruct, bea_geo,...
            bea_phys, bea_load);
    end
else if strcmp(rc_scenario,'inner')
        display('Outer race control detected. Re-running solver.');
        rc_scenario='outer';
        loadSolv=@(val)solvStitcher(obj, bea_geo, bea_phys, bea_load,...
            rc_scenario, val);
        [solVec,fval,exitflag,output]=fsolve(loadSolv,init,...
            options);
        resStruct = addSigmaMaxAlpha(obj, resStruct, bea_geo,...
            bea_phys, bea_load);
end
end

resStruct.raceControl=rc_scenario;

dispArray=[resStruct.delta_a, resStruct.delta_ry...
    resStruct.delta_rz resStruct.M_y...
    resStruct.M_z];
end
