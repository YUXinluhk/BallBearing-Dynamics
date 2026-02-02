function res = innerGeometryError(X_1, X_2, A_1,...
    A_2, delta_i, f_i, D)
res = (A_1-X_1).^2+(A_2-X_2).^2-((f_i-0.5)*D+delta_i).^2;
end
