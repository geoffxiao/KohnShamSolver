% Calculate the occupancy of each Kohn-Sham orbital

% Inputs:
%   Kohn-Sham orbital energies
%   Kohn-Sham orbitals
%   # of electrons
%   angular momentum values to search
%   radial axis

% Outputs
%   f_nl occupancy factor
%   electron density

function out = CalcOcc(E_allowed, wavefunctions_allowed, num_electrons, l_tot, r_axis)

    % Given allowed energies, let's find the occupancy of each allowed
    % energy!
    
    % Aufbau = occupy the lowest energy level first
    
    % Ignore spins = 2 electrons can occupy the same level
    
    % Total number of electrons = sum of all occupied states
    
    
    sorted_energy = E_allowed; % Will set to 0 < f_nl < 1 for wavefunction occupancy
    % Find the lowest energy wavefunction {n,l} order in terms of 1,2,3...
    temp = [];
    for i = 0 : l_tot
    	temp = [temp; E_allowed{i + 1}];    
    end
    [sorted, indx] = sort(temp);
    
    % Find the {n,l} wavefunctions that are occupied
    for c = 1 : numel(sorted)
        for l = 0 : l_tot % l number
            % n number for energy level
            n = find(E_allowed{l + 1} == sorted(c));
            if numel(n) > 0 % this l number has the lowest energy level
                sorted_energy{l + 1}(n) = c; % label that energy
            end
        end
    end

    % create f_nl occ factor
    f_nl = E_allowed;
    for l = 0 : l_tot
        for n = 1 : numel(f_nl{l + 1})
            f_nl{l + 1}(n) = 0;
        end
    end
    
    % fill the states from the lowest first
    tot = 0;
    for c = 1 : numel(sorted)
        for l = 0 : l_tot % l number
            n = find(sorted_energy{l + 1} == c); % find the current lowest energy level
            if numel(n) > 0 % found the current l number
                % find f_nl
                to_fill = (num_electrons/2) - tot; % number of electrons we need
                max_fill = 2 * l + 1;
                if to_fill > max_fill % need to fill all degen states
                    f_nl{l + 1}(n) = 1; % all states filled
                    tot = tot + max_fill;
                elseif to_fill > 0
                    f_nl{l + 1}(n) = to_fill / (2 * l + 1);
                    tot = tot + (2 * l + 1) * f_nl{l + 1}(n);
                end
            end    
        end
    end
    
    % loop through all energies to calc electron density
    ElectronDensity = 0;
    for l = 0 : l_tot
        for n = 1 : numel(f_nl{l + 1})
            ElectronDensity = ElectronDensity + ...
                2 * f_nl{l + 1}(n) * (2 * l + 1) .* ...
                (wavefunctions_allowed{l + 1}(:,n)).^2 ./ (4 * pi * r_axis.^2);
        end
    end
    
    out = cell(2,1);
    out{1} = ElectronDensity;
    out{2} = f_nl;
    
end