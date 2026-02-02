function obj = setGeometry(obj, D, d_m, alpha_free, r_i, r_o, z,...
    Delta_P_D, psi_0)
% setGeometry Sets the BallBearing object geometry used in
% displacement and stiffness calculations
obj.geometry.D=D;
obj.geometry.d_m=d_m;
obj.geometry.alpha_free=alpha_free;
obj.geometry.r_i=r_i;
obj.geometry.r_o=r_o;
obj.geometry.z=z;
obj.geometry.Delta_P_D=Delta_P_D;
obj.geometry.psi=genPsi(z, psi_0);
obj.geometry.f_i=r_i/D;
obj.geometry.f_o=r_o/D;
obj.geometry.R_i=calcRadLoc(d_m, obj.geometry.f_i, D,...
    alpha_free);
obj.geometry.gamma_tick=D/d_m; %p.308
obj.geometry.set=1;

obj.physical.set=0; %obj.physical needs updating after this
end
