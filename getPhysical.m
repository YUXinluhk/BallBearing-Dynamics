function phys = getPhysical(obj) %no pun intended
%generates a cell array of the input values for setPhysical
%so the properties can be updated by calling
%BallBearing.setPhysiscal(Ballbearing.getPhysical{:})
phys = cell(1,5);
phys{1} = obj.physical.rho_ball;
phys{2} = obj.physical.E_I;
phys{3} = obj.physical.E_II;
phys{4} = obj.physical.xi_I;
phys{5} = obj.physical.xi_II;
end
