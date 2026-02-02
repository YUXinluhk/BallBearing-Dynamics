function res=calcStiffness(obj, loads, init, h_force)
%for loads and init see calcDisplacement. h_force is the
%force step width for the numerical derivation.
[a,init]=calcDisplacement(obj, loads, 'outer', init);
res=ones(3,3);
if isempty(h_force)
    h_force=sqrt(eps);
end
for i=1:3;
    load_pre=loads;
    load_pre(i)=load_pre(i)-h_force/2;
    load_aft=loads;
    load_aft(i)=load_aft(i)+h_force/2;
    [~,~,pre_disp]=calcDisplacement(obj, load_pre, a.raceControl, init);
    [~,~,aft_disp]=calcDisplacement(obj, load_aft, a.raceControl, init);
    res(:,i)=h_force./(aft_disp(1:3)-pre_disp(1:3));
end
%res(abs(res)>1e17)=inf; %For better readability, stiffnesses
%greater than the calculation precision can be set to inf
end
