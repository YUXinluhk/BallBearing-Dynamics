function [Q, sigma_max, a, E, K, b] = calcStress(obj, delta, rho, phys)
%This function calculates everything about a Herzian ellipsoid
%that the class needs in various places. In order: force,
%peak stress, semi-major axis, elliptical integral, stiffness,
%and semi-minor axis.
persistent kappa
sum_rho = curvSum(rho);

%precise (see way down) determines the way that the kappa value
% of the ellipse is calculated:
% 0: By correlation as reccomended by Harris
% 1: Numerically
% Numerical calculation is currently very slow and generally
% not worth the extra calculation time
if precise
    if isempty(kappa)
        R_y = 1./(rho(:,1)+rho(:,3));
        R_x = 1./(rho(:,2)+rho(:,4));
        kappaInit = ellipseCorr(R_x, R_y);
    else kappaInit=kappa;
    end

    kapSolv=@(kap)hertzError(kap,rho);
    options = optimoptions('fsolve','Display','none', ...
        'MaxFunEvals', 1e9,...
        'MaxIter', 1e9, ...
        'Jacobian', 'on',...
        'FunValCheck','on',...
        'Algorithm','levenberg-marquardt',...
        'TolFun',1e-10,...
        'TolX', 1e-10 ...
        );
    % trust-region-reflective 'FinDiffRelStep', 1e-6,...
    % trust-region-dogleg

    [kappa,~,exitflag]=fsolve(kapSolv, kappaInit, options);
    if ~exitflag
        error('Geometry Error: Hertzian ellipsoid did not converge');
    end
    if ~isreal(kappa)
        error('Kappa not real')
    end
else
    %                 F_rho = curvDiff(rho);
    R_y = 1./(rho(:,1)+rho(:,3));
    R_x = 1./(rho(:,2)+rho(:,4));

    % Ensure Rx >= Ry for valid Hertzian formulas (kappa >= 1)
    % If Rx < Ry, swap input radii, calculate, then swap output axes a/b
    swapped = (R_x < R_y);
    if any(swapped)
        temp = R_x;
        R_x(swapped) = R_y(swapped);
        R_y(swapped) = temp(swapped);
    end

    kappa = ellipseCorr(R_x, R_y);
end
[F, E] = ellipke(1-1./kappa.^2); %Non-standard definition
%of the input parameter:
% m{matlab}=k^2{wikipedia}=1-1/kappa^2{Jones/Harris}

delta_star=2*F/pi.*(pi./(2*kappa.^2.*E)).^(1/3);
a_star=((2*kappa.^2.*E)/pi).^(1/3);
b_star=((2*E)./(pi*kappa)).^(1/3);

% Swap axes back if needed
if any(swapped)
    temp_a = a_star;
    a_star(swapped) = b_star(swapped);
    b_star(swapped) = temp_a(swapped);
end
%Gupta has 2 instead of pi in denominator in his book,
%but pi is correct
one_over_E_tick=elasticity(phys.xi_I, phys.xi_II,...
    phys.E_I, phys.E_II);
Q = engPower(delta./delta_star.*2./sum_rho,(3/2))...
    .*(2*sum_rho)./(3*one_over_E_tick);
a = a_star.*engPower(3/2*Q./sum_rho*one_over_E_tick,(1/3));
b = b_star.*engPower(3/2*Q./sum_rho*one_over_E_tick,(1/3));
sigma_max=(3*Q)./(2*pi*a.*b);
K=pi*kappa*2./(one_over_E_tick).*engPower(2*E./(9*sum_rho.*F.^3),0.5);
end
