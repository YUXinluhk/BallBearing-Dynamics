function F_rho = curvDiff(rho)
F_rho = ((rho(:,1)-rho(:,2))+(rho(:,3)-rho(:,4)))./sum(rho,2);
end
