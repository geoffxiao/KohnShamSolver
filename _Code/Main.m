% Kohn-Sham solver for single atoms
% Geoffrey Xiao, Cooper Voigt, Evan Quain
% MATE 460/I699 Final Project

% To run... Inputs:
%   Only need to change parameters in Setup
%   STRING = what to save the file as
%   Z_num = # of protons
%   num_electrons = # of electrons
%   r_axis = the radial axis
%   f_factor = mixing factor
%   N_basis_func = # of basis functions to use
%   l_tot = # of l values (orbital shapes s, p, d, f) to search
%   convergence_criterion = when to stop self-consistency


% Outputs:
%   E_allowed = the energies of the Kohn-Sham orbitals
%   f_nl = occupancy of the orbitals
%   E_total = total energy of the system, can compare with ionization energies in NIST database
%   ElectronDensity = electron density (m^-3)
%   RadialElectronDensity = electron density at a radius r (m^-1)


clear all;
clc;

%% Setup

STRING = sprintf('H'); % name to save 
r_axis = [1e-10 : 0.001 : 5]'; % radial axis, can't start at 0 b/c diverging

% What type of atom?
Z_num = 2; % # of protons
num_electrons = Z_num; % # of electrons

% Convergence variables
f_factor = 0.25; % mixing function for self-consistent
N_basis_func = 200; % num of basis functions to use

l_tot = 3; % number of orbital types (s,p,d,f) to calc... angular momentum searching
convergence_criterion = 1e-4; % 1e-4 is more than fine


%%----------------------------------------------------------------------------------------------------------------------%%
%%----------------------------------------------------------------------------------------------------------------------%%
%% DO NOT MODIFY BELOW THIS LINE
%%----------------------------------------------------------------------------------------------------------------------%%
%%----------------------------------------------------------------------------------------------------------------------%%




%% Set up Basis Set Expansion
Lr = r_axis(end); % for basis function formation
% Set up the basis functions, use the particle in a box basis set (sine
% functions)
BasisFunctions = zeros(numel(r_axis), N_basis_func); 
for i = 1 : N_basis_func
    BasisFunctions(:,i) = sqrt(2/Lr) * sin(i*pi*(r_axis-r_axis(1))/Lr);
end

% Secular_Matrix KE Part
SecularMatrix_T = zeros(numel(r_axis), N_basis_func, l_tot + 1);
for l = 0 : l_tot
    for i = 1 : N_basis_func
        SecularMatrix_T(:,i,l + 1) = ...
            sqrt(1/(2*Lr)) .* (i*pi/Lr)^2 .* sin(i*pi*r_axis/Lr) + ...
            ((l*(l+1))./(2*r_axis.^2)) .* sqrt(2/Lr) .* sin(i*pi*r_axis/Lr);
    end
end

%% START self-consistency
iteration_counter = 1;

% first guesses
Error = inf;
ElectronDensity_curr = ( num_electrons^4 / (64 * pi) ) * exp( -num_electrons * r_axis / 2); % current guess
Errors = [];
E_calc = [];

V_nuc = (-Z_num ./ r_axis); % nuclear potential function

while( Error > convergence_criterion )
 
    [E_exchange, V_exchange, E_correlation, V_correlation] = XCEnergy( ElectronDensity_curr ); % calc exchange corr energies
    V_hartree = HartreeEnergy( ElectronDensity_curr, r_axis ); % calc hartree energy
    V_curr = V_exchange + V_correlation + V_hartree + V_nuc; % the total potential energy
    V_curr = V_curr - V_curr(end);
    
    shro_sol = SolveShrodinger( BasisFunctions, SecularMatrix_T, V_curr, l_tot, r_axis ); % solve shrodinger equation
    E_allowed = shro_sol{1}; % allowed energies
    wavefunctions_allowed = shro_sol{2}; % The radial part * r
    
    occ = CalcOcc(E_allowed, wavefunctions_allowed, num_electrons, l_tot, r_axis); % calc occupancy of each orbital
    ElectronDensity_new = occ{1}; % electron density, found a new electron density
    f_nl = occ{2}; % occ of each wavefunction

    % energy of this electron density
    [E_total, E_k, E_nuc, E_hartree, E_exchange, E_correlation] ...
    = TotalEnergy( ElectronDensity_new, r_axis, E_exchange, V_exchange, E_correlation, ...
    V_correlation, V_hartree, V_nuc, f_nl, E_allowed, l_tot);
    
    iteration_counter = iteration_counter + 1;
    
    % next electron density guess is mix of the old and the new
    ElectronDensity_curr = ElectronDensity_new * f_factor + ElectronDensity_curr * (1 - f_factor); % new ElectronDensity
    Error = sum(abs( ElectronDensity_curr - ElectronDensity_new )) % has it converged?
    
    % The errors
    Errors = [Errors; Error];
    E_calc = [E_calc; E_total];
    
end

%% Calc electron densities
ElectronDensity = ElectronDensity_curr;
% find the radial electron probability density func
RadialElectronDensity = 0;
for l = 0 : l_tot
    for n = 1 : numel(f_nl{l + 1})
        RadialElectronDensity = RadialElectronDensity + ...
            2 * f_nl{l + 1}(n) * (2 * l + 1) .* ...
            (wavefunctions_allowed{l + 1}(:,n)).^2;
    end
end

% the various energies
% E_total, E_k, E_nuc, E_hartree, E_exchange, E_correlation;

%% save data
save(sprintf('%s.mat',STRING),'E_allowed','f_nl','RadialElectronDensity',...
'wavefunctions_allowed','E_total','r_axis','N_basis_func','ElectronDensity')