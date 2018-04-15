% Solve Poisson's equation to calculate Hartree potential energy

% Inputs:
%   Electron density
%   Radial axis

% Outputs:
%   Hartree potential energy


function V_hartree = HartreeEnergy( ElectronDensity, r_axis )

    h_spacing = r_axis(2) - r_axis(1); % radial axis spacing
    Lr = r_axis(end);

    % Poisson Equation
    s = -4 * pi * r_axis .* ElectronDensity;
    % U_e = cumtrapz( r_axis, cumtrapz(r_axis,s) );
    U_e = zeros(numel(r_axis),1);
    U_e(1) = 0;
    U_e(2) = h_spacing/200;
    for i = 2 : numel(U_e) - 1
        U_e(i + 1) = 2 * U_e(i) - U_e(i - 1) + (h_spacing^2 / 12)*(s(i + 1) + 10 * s(i) + s(i - 1));
    end
    
    TotalCharge = 4 * pi * trapz(r_axis, r_axis.^2 .* ElectronDensity);
    U_e = U_e + r_axis .* (TotalCharge - U_e(end)) ./ Lr;
    V_hartree = U_e ./ r_axis;

end