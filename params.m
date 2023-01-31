
function [alpha, nu, delta, beta, qss, exit, psi, zeta, tax, opp_cost, zss, knum, kgrid, k_enter, bknum, bkgrid, bk_enter, enum, epsgrid, pie] = params()


% standard parameters
% K share        % N share     % cap dep         % HH disc        % debt         % exit          % leisure
alpha = 0.28;    nu = 0.60;    delta = 0.069;    beta  = 0.96;    qss = beta;    exit  = 0.1;    psi   = 2.11;

% fun parameters
zeta = .91;     tax = 0.05;    opp_cost = 0.0;   zss = 1.0;   % 91/68   %05/066/1

% capital choice grid
klow = 0.01; khigh = 7; knum = 75;
kgrid = linspace(klow, khigh, knum)';
k_enter = .125;

% bond choice grid
bklow = -15.0; bkhigh = 1.5; bknum = 75;                    % K/O = 2.28    .122      .337     .372
bkgrid = linspace(bklow, bkhigh, bknum)';
bk_enter = 0; %min(abs(bkgrid));

% idio shock structure
%mean          %N_e        % exp autoreg   %std innov
meane = 0.0;   enum = 7;    rho = 0.665;   stdinnov = 0.057;
[epsgrid, pie] = tauchen2(meane, stdinnov, rho, 2.25, enum);
epsgrid = exp(epsgrid);

pi_zeta = [0.9765 0.0235; 0.3125 0.6875];