%backwards
% needs:      T, zeta, bail_rate, wage, ldf(q), to work    

function [y_val, n, kft, bkft, w_vec, q_vec, bailout_t] = backwards(T, zeta_vec, bail_rate, w_vec, q_vec);

[alpha, nu, delta, beta, qss, exit, psi, zeta, tax, opp_cost, zss, knum, kgrid, k_enter, bknum, bkgrid, bk_enter, enum, epsgrid, pie] = params();

% dec rules are t dependent 
kft = zeros(knum, bknum, enum, T);
bft = kft;
bkft = kft;
bailout_t = kft;

% things are not t dependent 
n = zeros(knum, bknum, enum);
y_val = n;
cash = n;
dive = n;
invest = n;


for t = 1:T

for eps_i = 1:enum
    for bk_i = 1:bknum
        for k_i = 1:knum

            bval = bkgrid(bk_i)*kgrid(k_i);

            n(k_i, bk_i, eps_i) = ( (nu*zss*epsgrid(eps_i)*kgrid(k_i)^(alpha)) / ((1+tax)*w_vec(1,t)) ) ^ (1 / (1-nu));

            y_val(k_i, bk_i, eps_i) = zss*epsgrid(eps_i)*kgrid(k_i)^(alpha)*n(k_i, bk_i, eps_i)^(nu);

            cash(k_i, bk_i, eps_i) = y_val(k_i, bk_i, eps_i) - opp_cost - (1+tax)*w_vec(1,t)*n(k_i, bk_i, eps_i) + (1-delta)*kgrid(k_i) - bval; % (1-bail)*bval


            % FOC of V wrt k'
            for j_eps = 1:enum
                kft(k_i, bk_i, eps_i, t) = kft(k_i, bk_i, eps_i, t) + pie(eps_i,j_eps)*alpha/(1-nu)*((zss*nu/(1+tax)*w_vec(1,t))^(nu/(1-nu))*epsgrid(j_eps)^(1/(1-nu))-(1+tax)*w_vec(1,t)*(zss*nu*epsgrid(j_eps)/(1+tax)*w_vec(1,t))^(1/(1-nu)));
            end
            kft(k_i, bk_i, eps_i, t) = ((1/beta-1+delta)/(kft(k_i, bk_i, eps_i, t)))^((1-nu)/(alpha+nu-1));

            % zero dive
            bft(k_i, bk_i, eps_i, t) = (1/q_vec(1,t)) * (kft(k_i, bk_i, eps_i, t) - cash(k_i, bk_i, eps_i));



            % bailout conditions
            if  bval > 0 && bft(k_i, bk_i, eps_i, t) > 0  && eps_i > 0.5
                bailout_t(k_i, bk_i, eps_i, t) = bail_rate*bval;
                cash(k_i, bk_i, eps_i) = cash(k_i, bk_i, eps_i) + bailout_t(k_i, bk_i, eps_i, t);
                bft(k_i, bk_i, eps_i, t) = bft(k_i, bk_i, eps_i, t) - bailout_t(k_i, bk_i, eps_i, t);
            end


            % col const
            if bft(k_i, bk_i, eps_i, t) > kgrid(k_i)*zeta_vec(1,t)
                bft(k_i, bk_i, eps_i, t) = kgrid(k_i)*zeta_vec(1,t);
                kft(k_i, bk_i, eps_i, t) =  cash(k_i, bk_i, eps_i) + q_vec(1,t)*bft(k_i, bk_i, eps_i, t);
            end

            bkft(k_i, bk_i, eps_i, t) = bft(k_i, bk_i, eps_i, t) / kft(k_i, bk_i, eps_i, t); % leverage tomorrow
        end
    end
end
end