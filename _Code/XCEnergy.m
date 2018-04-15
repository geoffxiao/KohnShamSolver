% Calculate exchange correlation energies

% Inputs:
%   Electron density

% Outputs:
%   energy densities
%   energy potential energies


function [E_x, V_x, E_c, V_c] = XCEnergy( ElectronDensity )

    % XC energies    
    c = -(3/(2*pi))^(2/3);
    rs = ( 3 ./ (4*pi*ElectronDensity) ).^(1/3);
    
    E_x = (3/4) * (c./rs);
    V_x = c ./ rs;
    
    a = 0.0621814; b = 3.72744; c = 12.9352; x0 = -0.10498;
    x = sqrt(rs);
    q = sqrt(4*c - b^2);
    
    P = @(x) x.^2 + b * x + c;
    D = @(x) rs ./ P(x);
    Y = @(x) atan(q ./ (2*x+b));
    W = @(x) (x - x0).^2 ./ P(x);
    
    f1 = 2*b/q;
    f2 = -b*x0/P(x0);
    f3 = 2*(b + 2*x0)*f2 / q;
    
    E_c = (a/2) * ( log(D(x)) + f1 * Y(x) + f2 * log(W(x)) + f3 * Y(x) );
    E_c(1) = 0;
    V_c = E_c - (a/6) * ( c*(x-x0) - b*x*x0 ) ./ ( (x - x0) .* P(x) );
    V_c(1) = 0;
       
end