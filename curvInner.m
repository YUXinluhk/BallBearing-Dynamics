function rho_ci = curvInner(cos_alp_i, f_i, D, d_m)
%calculates the curvature at the inner raceway contact based on
%Jones' derivation.
gamma = D*cos_alp_i/d_m;
%I: Ball
%II: Raceway
%1: Ellipsoid axis in raceway curvature direction (typically
%minor axis)
%2: Ellipsoid axis in grove curvature direction (typically
%major axis)
rho_I1 = 2/D+0*gamma;
rho_I2 = 2/D+0*gamma;
rho_II1 = (2/D)*(gamma./(1-gamma));
rho_II2 = -1/(f_i*D)+0*gamma;
rho_ci=[rho_I1 rho_I2 rho_II1 rho_II2];
end
