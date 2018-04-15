% Calculate the total energy

% Inputs:
%   Energy densities 
%   Electron density

% Outputs:
%   the kinetic, nuclear, hartree, exchange, and correlation energies
%   the total energy of the system


function [E_total, E_k, E_nuc, E_hartree, E_exchange, E_correlation] ...
    = TotalEnergy( ElectronDensity, r_axis, E_exchange, V_exchange, E_correlation, V_correlation,...
    V_hartree, V_nuc, f_nl, E_allowed, l_tot)

    V_curr = V_exchange + V_correlation + V_hartree + V_nuc;
    
    V_tilde = (E_exchange - V_exchange) + (E_correlation - V_correlation) - 0.5 * V_hartree;
    g = V_tilde .* ElectronDensity;
    E_sum = 0;
    for l = 0 : l_tot
        for n = 1 : numel(f_nl{l + 1})
            E_sum = E_sum + ...
                2 * f_nl{l + 1}(n) * (2 * l + 1) .* ...
                E_allowed{l + 1}(n);
        end
    end

    E_total = E_sum + trapz(r_axis, 4 * pi * r_axis.^2 .* g);
    E_k = E_sum - trapz(r_axis, 4 * pi * r_axis.^2 .* V_curr .* ElectronDensity);
    E_nuc = trapz(r_axis, 4 * pi * r_axis.^2 .* V_nuc .* ElectronDensity);
    E_hartree = 0.5 * trapz(r_axis, 4 * pi * r_axis.^2 .* V_hartree .* ElectronDensity);
    E_exchange = trapz(r_axis, 4 * pi * r_axis.^2 .* E_exchange .* ElectronDensity);
    E_correlation = trapz(r_axis, 4 * pi * r_axis.^2 .* E_correlation .* ElectronDensity);

end