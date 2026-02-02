function A_2 = raceCenterDistRad(delta_ry, delta_rz, geo)
%Jones includes a radial clearance here, but Harris doesn't.
%Can be added without a problem.
A_2 = (geo.f_i+geo.f_o-1)*geo.D*cosd(geo.alpha_free) + ...
    delta_rz*cosd(geo.psi)+delta_ry*sind(geo.psi);
if ~isreal(A_2)
    disp('A_2 complex')
end
end
