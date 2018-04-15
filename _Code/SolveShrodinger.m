% Solve Shrodinger's equation

% Inputs:
%   basis functions
%   the secular matrix with the kinetic term
%   potential energy
%   l values to search
%   radial axis

% Outputs:
%   energies of the Kohn-Sham orbitals
%   Kohn-Sham orbitals


function out = SolveShrodinger( BasisFunctions, SecularMatrix_T, V_curr, l_tot, r_axis )

    % Solve Shrodinger's equation, get allowed energies and wavefunctions
    [~, N_basis_func] = size(BasisFunctions);
    SecularMatrix = zeros(N_basis_func, N_basis_func, l_tot);
    wavefunctions_allowed = cell(1, l_tot);
    E_allowed = cell(1, l_tot);
    
    % Basis set expansion method
    
    for l_num = 0 : l_tot % calc for diff. l numbers

        % Solve Shrodinger
        SecularMatrix_l_num = zeros(N_basis_func, N_basis_func);
        for i = 1 : N_basis_func % column
            for j = i : N_basis_func % row

                % find the hamiltonian operating on basis function i
                SecularMatrix_V = V_curr .* BasisFunctions(:,i);
                Hamiltonian_i = SecularMatrix_T(:,i,l_num + 1) + SecularMatrix_V;

                % store in the secular matrix
                SecularMatrix_l_num(j, i) = trapz(r_axis, BasisFunctions(:,j) .* Hamiltonian_i);
                SecularMatrix_l_num(i, j) = SecularMatrix_l_num(j, i);

            end
        end

        SecularMatrix(:,:,l_num + 1) = SecularMatrix_l_num;

        % coefficients to make the wavefunctions
        % eigenenergies
        [coeffs_l_num, energies_l_num] = eig(SecularMatrix_l_num);
        wavefunctions_allowed_l_num = BasisFunctions * coeffs_l_num;
        E_allowed_l_num = diag(energies_l_num);
        E_allowed_l_num = E_allowed_l_num(E_allowed_l_num < max(V_curr)); % E must be less than max V -> Bound States
        wavefunctions_allowed_l_num = wavefunctions_allowed_l_num(:,1 : numel(E_allowed_l_num));

        wavefunctions_allowed{l_num + 1} = wavefunctions_allowed_l_num;
        E_allowed{l_num + 1} = E_allowed_l_num;

    end
    
    out = cell(2,1);
    out{1} = E_allowed;
    out{2} = wavefunctions_allowed;
    
end