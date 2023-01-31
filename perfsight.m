clear
close all
load pfinputs

[alpha, nu, delta, beta, qss, exit, psi, zeta, tax, opp_cost, zss, knum, kgrid, k_enter, bknum, bkgrid, bk_enter, enum, epsgrid, pie] = params();

zeta = 0.91;

T=2;
zeta_vec = ones(1,T)*zeta;          % 91/68   %05/066/1
bail_rate = 0.0;
pf_tol = 1e-4;
pf_dist = 1;

c_guess = ones(1,T+1)*css;
c_estim = c_guess;
lambda = .95;

pfIter = 0;
while pf_dist > pf_tol
pfIter = pfIter + 1;


for t = 1:T
    w_vec(1,t) = psi*c_guess(t);  % psi * c
    q_vec(1,t) = beta * c_guess(t)/c_guess(t+1);      % beta * c/c'
end


% BACKWARDS needs: T, zeta, bail_rate, wage, ldf(q), to work   
[y_val, n, kft, bkft, w_vec, q_vec, bailout_t] = backwards(T, zeta_vec, bail_rate, w_vec, q_vec);


% FORWARDS needs: T, mu, yval, n, kf, bkf to work
[agg_k, agg_b, agg_n, agg_y, Gov_Et, w_vec, q_vec, dist] = forwards(T, mu, y_val, n, kft, bkft, w_vec, q_vec, bailout_t);

for t = 1:T
    c_estim(t) = agg_y(t) - agg_k(t+1) + (1-delta)*agg_k(t) - Gov_Et(t);
end

pf_dist = max(abs(c_guess(1:T) - c_estim(1:T)));

c_guess = c_guess * (lambda) + c_estim * (1-lambda);

end