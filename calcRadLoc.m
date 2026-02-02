function R_i = calcRadLoc(d_m, f_i, D ,alpha_free)
R_i = 0.5*d_m+(f_i-0.5)*D*cosd(alpha_free);
end
