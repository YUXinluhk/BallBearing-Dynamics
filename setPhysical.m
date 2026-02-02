function obj = setPhysical(obj, rho_ball, E_I, E_II, xi_I, xi_II)
if obj.geometry.set==0
    error(['Bearing geometry not set. Use method'...
        ' setGeometry before calling this method']);
else
    obj.physical.rho_ball=rho_ball;
    obj.physical.m_ball=...
        ballMass(obj.geometry.D,obj.physical.rho_ball);
    obj.physical.J_ball=ballInertiaMoment(obj.geometry.D,obj.physical.m_ball);
    obj.physical.E_I=E_I;
    obj.physical.E_II=E_II;
    obj.physical.xi_I=xi_I;
    obj.physical.xi_II=xi_II;
    obj.physical.set=1;
end
end
