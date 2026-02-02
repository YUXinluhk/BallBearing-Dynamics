function psi=genPsi(z,psi_0)
%devide a 360 deg circle into z equal parts starting at psi_0
psi=((0:(360/z):360-(360/z))+psi_0)';
end
