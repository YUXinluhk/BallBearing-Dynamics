function res = outerGeometryError(X_1, X_2, delta_o, f_o, D)
res = X_1.^2+X_2.^2-((f_o-0.5)*D+delta_o).^2;
end
