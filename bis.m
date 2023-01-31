
           function [agg_B, agg_C, agg_I, agg_K, agg_N, agg_Y, wss, mu] = bis()

% clear
% close all

% standard parameters
% K share        % N share     % cap dep         % HH disc        % debt         % exit          % leisure
alpha = 0.28;    nu = 0.60;    delta = 0.069;    beta  = 0.96;    qss = beta;    exit  = 0.1;    psi   = 2.11;

% fun parameters
zeta = .91;     tax = 0.05;    bail_rate = 0.0;   opp_cost = 0.0;   zss = 1.0;   % 91/68   %05/066/1

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialization & Bisection for pval = DU(c)
bisTol=1e-7; 
bisDistance=2*bisTol;

v0 = zeros(knum,bknum,enum);
bisIter = 0;
plow = .5; phigh = 4;

pval = plow;
[cash, dive, agg_C, agg_B, agg_K, agg_Y, agg_I, agg_N, Gov_E, total_bail, bailout, wss, v, mu, invest, kf, bf, check, agg_Cash, agg_Dive, bkf] = firm_dec(v0, pval, alpha, nu, delta, beta, qss, exit, psi, zeta, tax, bail_rate, zss, knum, kgrid, k_enter, bknum, bkgrid, bk_enter, enum, epsgrid, pie, opp_cost);
flow = pval - (1/agg_C);

s = sprintf( '  plow    p = %8.6f   cagg  = %8.4f   f_check = %8.6f', pval, agg_C, flow);
disp (s)

v0 = v;

pval = phigh;
[cash, dive, agg_C, agg_B, agg_K, agg_Y, agg_I, agg_N, Gov_E, total_bail, bailout, wss, v, mu, invest, kf, bf, check, agg_Cash, agg_Dive, bkf] = firm_dec(v0, pval, alpha, nu, delta, beta, qss, exit, psi, zeta, tax, bail_rate, zss, knum, kgrid, k_enter, bknum, bkgrid, bk_enter, enum, epsgrid, pie, opp_cost);
fhigh = pval - (1/agg_C);

s = sprintf( ' phigh     p = %8.6f   cagg  = %8.4f   f_check = %8.6f ', pval, agg_C,  fhigh);
disp (s)

if (flow*fhigh > 0.0)
    s = sprintf ( ' could not bisect between (p0, p1) = (%6.3f, %6.3f)   (pval0, pval1) = (%8.4f, %8.4f) ', plow, phigh, flow, fhigh);
    disp(s)
else
    
    while (bisDistance > bisTol)

        pval = (plow + phigh)/2.0;
        [cash, dive, agg_C, agg_B, agg_K, agg_Y, agg_I, agg_N, Gov_E, total_bail, bailout, wss, v, mu, invest, kf, bf, check, agg_Cash, agg_Dive, bkf] = firm_dec(v0, pval, alpha, nu, delta, beta, qss, exit, psi, zeta, tax, bail_rate, zss, knum, kgrid, k_enter, bknum, bkgrid, bk_enter, enum, epsgrid, pie, opp_cost);
        f = pval - (1/agg_C);

        if (f*flow > 0.0)
            plow = pval; flow = f;
        else
            phigh = pval; fhigh = f;
        end

        bisDistance  = abs(phigh - plow);
        bisIter = bisIter + 1;

        disp ( ' ' )
        s = sprintf( ' %4d     p = %8.6f   wage = %8.6f   cagg  = %8.4f ', bisIter, pval, wss, agg_C);
        disp(s)

        v0 = v;
        mu0 = mu;

    end
end

save ssdist




% %%%%%%%%%%% K X EPS %%%%%%%%%%%%%%%%%
% KEmu = zeros(knum, enum);
% [Emat, Kmat] = meshgrid(epsgrid, kgrid);
% 
% for k_i = 1:knum
%     for eps_i = 1:enum
%         KEmu(k_i, eps_i)= sum(mu(k_i,:,eps_i));
%     end
% end
% 
% figure
% 
% mesh(Kmat, Emat, KEmu)
% set(gca,'FontSize',15);
% set(get(gca,'Xlabel'),'FontSize',10)
% set(get(gca,'Ylabel'),'FontSize',10)
% hold on
% xlabel( ' Capital ', 'Fontsize', 14 )
% ylabel( ' Idio. productivity ','Fontsize', 14  )
% title ( ' Capital-productivity distribution ','Fontsize', 14  )
% 
% 
% %%%%%%%%%%% SS Dist %%%%%%%%%%%%%%%%%
% medmu = mu(:,:,4);
% [Kmat, BKmat] = meshgrid (kgrid, bkgrid);
% 
% figure
% 
% mesh(Kmat, BKmat, medmu')
% set(gca,'FontSize',30);
% set(get(gca,'Xlabel'),'FontSize',20)
% set(get(gca,'Ylabel'),'FontSize',20)
% hold on
% xlabel( ' Capital ', 'Fontsize', 14 )
% ylabel( ' Leverage ','Fontsize', 14  )
% title ( ' Steady-state distribution: median productivity ','Fontsize', 14  )
% 
% 

% %%%%%%%%%%% Value %%%%%%%%%%%%%%%%%
% [Kmat, BKmat] = meshgrid (kgrid, bkgrid);
% medv = v(:,:,4);
% 
% figure
% 
% mesh(Kmat, BKmat, medv')
% set(gca,'FontSize',30);
% set(get(gca,'Xlabel'),'FontSize',20)
% set(get(gca,'Ylabel'),'FontSize',20)
% hold on
% xlabel( ' Capital ', 'Fontsize', 14 )
% ylabel( ' Leverage ','Fontsize', 14  )
% title ( ' Value function ','Fontsize', 14  )