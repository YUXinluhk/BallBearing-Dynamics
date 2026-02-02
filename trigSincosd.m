function [cos_alp_i, cos_alp_o, sin_alp_i, sin_alp_o]=...
    trigSincosd(obj, A_1, A_2, X_1, X_2, delta_i, delta_o, D,...
    f_i, f_o)
cos_alp_i = cosICAngle(A_2, X_2, f_i, D, delta_i);
cos_alp_o = cosOCAngle(X_2, f_o, D, delta_o);
sin_alp_i = sinICAngle(A_1, X_1, f_i, D, delta_i);
sin_alp_o = sinOCAngle(X_1, f_o, D, delta_o);
end
