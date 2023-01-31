function [cash, dive, agg_C, agg_B, agg_K, agg_Y, agg_I, agg_N, Gov_E, total_bail, bailout, wss, v, mu, invest, kf, bf, check, agg_Cash, agg_Dive, bkf] = firm_dec(v0, pval, alpha, nu, delta, beta, qss, exit, psi, zeta, tax, bail_rate, zss, knum, kgrid, k_enter, bknum, bkgrid, bk_enter, enum, epsgrid, pie, opp_cost);


wss = psi/pval;
v = zeros(knum, bknum, enum);
ev = v;
tv = v;
vTol = 1e-7;
vDist = 2*vTol;

kf = zeros(knum, bknum, enum);
bf = kf;
bkf = kf;
n = kf;
y_val = kf;
cash = kf;
dive = kf;
invest = kf;
bailout = kf;

check = kf;


for eps_i = 1:enum
    for bk_i = 1:bknum
        bkval = bkgrid(bk_i);
        for k_i = 1:knum
            kval = kgrid(k_i);
            bval = bkval*kval;

            n(k_i, bk_i, eps_i) = ( (nu*zss*epsgrid(eps_i)*kval^(alpha)) / ((1+tax)*wss) ) ^ (1 / (1-nu));

            y_val(k_i, bk_i, eps_i) = zss*epsgrid(eps_i)*kval^(alpha)*n(k_i, bk_i, eps_i)^(nu);

            cash(k_i, bk_i, eps_i) = y_val(k_i, bk_i, eps_i) - opp_cost - (1+tax)*wss*n(k_i, bk_i, eps_i) + (1-delta)*kval - bval; % (1-bail)*bval


            % FOC of V wrt k'
            for j_eps = 1:enum
                kf(k_i, bk_i, eps_i) = kf(k_i, bk_i, eps_i) + pie(eps_i,j_eps)*alpha/(1-nu)*((zss*nu/(1+tax)*wss)^(nu/(1-nu))*epsgrid(j_eps)^(1/(1-nu))-(1+tax)*wss*(zss*nu*epsgrid(j_eps)/(1+tax)*wss)^(1/(1-nu)));
            end

            kf(k_i, bk_i, eps_i) = ((1/beta-1+delta)/(kf(k_i, bk_i, eps_i)))^((1-nu)/(alpha+nu-1));

            bf(k_i, bk_i, eps_i) = (1/qss) * (kf(k_i, bk_i, eps_i) - cash(k_i, bk_i, eps_i));



            % bailouts
            if  bval > 0 && bf(k_i, bk_i, eps_i) > 0  && eps_i > 6.5
                bailout(k_i, bk_i, eps_i) = bail_rate*bval;
                cash(k_i, bk_i, eps_i) = cash(k_i, bk_i, eps_i) + bailout(k_i, bk_i, eps_i);
                bf(k_i, bk_i, eps_i) = bf(k_i, bk_i, eps_i) - bailout(k_i, bk_i, eps_i);
            end


            % col const
            if bf(k_i, bk_i, eps_i) > kval*zeta
                check(k_i, bk_i, eps_i) = 1; % if check = 1, firm hit col contst

                bf(k_i, bk_i, eps_i) = kval*zeta;

                kf(k_i, bk_i, eps_i) =  cash(k_i, bk_i, eps_i) + qss*bf(k_i, bk_i, eps_i);
            end

            invest(k_i, bk_i, eps_i) =  kf(k_i, bk_i, eps_i) - kgrid(k_i);
            bkf(k_i, bk_i, eps_i) = bf(k_i, bk_i, eps_i) / kf(k_i, bk_i, eps_i); % leverage tomorrow
            dive(k_i, bk_i, eps_i) = cash(k_i, bk_i, eps_i) - kf(k_i, bk_i, eps_i) + qss*bf(k_i, bk_i, eps_i);
        end
    end
end

kmin = min(min(min(kf))); kmax = max(max(max(kf)));
bkmin = min(min(min(bkf))); bkmax = max(max(max(bkf)));
bfmin = min(min(min(bf))); bfmax = max(max(max(bf)));
cashmin = min(min(min(cash))); cashmax = max(max(max(cash)));
Dmin = min(min(min(dive))); Dmax = max(max(max(dive)));

print0 = sprintf ( '\n\t(kmin, kmax) = (%2.8f, %2.8f), (bfmin, bfmax) = (%2.8f, %2.8f), (bkmin, bkmax) = (%2.8f, %2.8f), \n\t(cashmin, cashmax) = (%2.8f, %2.8f), (Dmin, Dmax) = (%2.8f, %2.8f)', ...
    kmin, kmax, bfmin, bfmax, bkmin, bkmax, cashmin, cashmax, Dmin, Dmax);
disp(print0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vIter = 0;
while (vDist > vTol)
    vIter = vIter + 1;

    %%%%%%%%%%%%%%%%%%
    %   set up EV    %
    %%%%%%%%%%%%%%%%%%

    for bf_i = 1:1:bknum
        for kf_i = 1:1:knum
            for eps_i = 1:1:enum
                pie_vec = pie(eps_i, :);
                temp_ev = 0.0;
                for eps_j = 1:1:enum
                    temp_ev = temp_ev + pie_vec(eps_j)*v(kf_i, bf_i, eps_j);
                end
                ev(kf_i, bf_i, eps_i) = temp_ev;

            end
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   THE value fn iter    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    for eps_i = 1:enum
        evbk = ev(:, :, eps_i);
        for bf_i = 1:bknum
            for kf_i = 1:knum

                bf_val = bf(kf_i, bf_i, eps_i);
                kf_val = kf(kf_i, bf_i, eps_i);
                bk_val = bf_val/kf_val;

                exit_v = exit*cash(kf_i, bf_i, eps_i);

                [k_ind, k_weight] = linear(knum, kgrid, kf_val);
                [bk_ind, bk_weightw] = linear(bknum, bkgrid, bk_val);

                cont_v = bk_weightw*evbk(:, bk_ind) + (1-bk_weightw)*evbk(:, bk_ind+1);
                cont_v = k_weight*cont_v(k_ind) + (1-k_weight)*cont_v(k_ind+1);
                cont_v = (1-exit)*beta*cont_v;

                tv(kf_i, bf_i, eps_i) = exit_v + cont_v;
            end
        end
    end


    vDist = max(max(max(abs(v - tv))));
    v = tv;


    kmin = min(min(min(kf))); kmax = max(max(max(kf)));
    bkmin = min(min(min(bkf))); bkmax = max(max(max(bkf)));
    cashmin = min(min(min(cash))); cashmax = max(max(max(cash)));
    Dmin = min(min(min(dive))); Dmax = max(max(max(dive)));

    if (mod(vIter,10) == 0)
        print0 = sprintf ( 'At (q, w) = (%2.4f, %2.4f), \n\titeration %4d, ||tv-v|| = %2.8f, \n\t(kmin, kmax) = (%2.8f, %2.8f), (bkmin, bkmax) = (%2.8f, %2.8f), \n\t(cashmin, cashmax) = (%2.8f, %2.8f), (Dmin, Dmax) = (%2.8f, %2.8f)', ...
            qss, wss, vIter, vDist, kmin, kmax, bkmin, bkmax, cashmin, cashmax, Dmin, Dmax);
        disp(print0);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                          EXITS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% firm distribution

distTol = 1e-11;
dist_error = 1;

%Initial guess
k_0 = 1;
bk_0 = 0;

stat_dist_eps = pie^1000;
stat_dist_eps = stat_dist_eps(1, :);


mu = zeros(knum,bknum,enum);
[k_ind,k_weight] = linear(knum, kgrid, k_0);
[bk_ind,bk_weight] = linear(bknum, bkgrid, bk_0);

for eps_i = 1:enum
    mu(k_ind,bk_ind,eps_i) = k_weight*bk_weight*stat_dist_eps(eps_i);
    mu(k_ind+1,bk_ind,eps_i) = (1-k_weight)*bk_weight*stat_dist_eps(eps_i);
    mu(k_ind,bk_ind+1,eps_i) = k_weight*(1-bk_weight)*stat_dist_eps(eps_i);
    mu(k_ind+1,bk_ind+1,eps_i) = (1-k_weight)*(1-bk_weight)*stat_dist_eps(eps_i);
end




while dist_error>distTol
    mu1 = zeros(knum, bknum, enum);
    for k_i = 1:knum
        for bk_i = 1:bknum
            for eps_i = 1:enum
                if mu(k_i,bk_i,eps_i) > 0
                    kf_val = kf(k_i,bk_i,eps_i);
                    bkf_val = bkf(k_i,bk_i,eps_i);

                    [k_ind,k_weight] = linear(knum, kgrid, kf_val);
                    [bk_ind,bk_weight] = linear(bknum, bkgrid, bkf_val);

                    % entrants
                    [k_ent_ind,k_ent_weight] = linear(knum, kgrid, k_enter);
                    [bk_ent_ind,bk_ent_weight] = linear(bknum, bkgrid, bk_enter);

                    for eps_j = 1:enum
                        mu1(k_ind,bk_ind,eps_j) = mu1(k_ind,bk_ind,eps_j) + k_weight*bk_weight*pie(eps_i,eps_j)*mu(k_i,bk_i,eps_i) *(1-exit);
                        mu1(k_ind+1,bk_ind,eps_j) = mu1(k_ind+1,bk_ind,eps_j) + (1-k_weight)*bk_weight*pie(eps_i,eps_j)*mu(k_i,bk_i,eps_i)*(1-exit);
                        mu1(k_ind,bk_ind+1,eps_j) = mu1(k_ind,bk_ind+1,eps_j) + k_weight*(1-bk_weight)*pie(eps_i,eps_j)*mu(k_i,bk_i,eps_i)*(1-exit);
                        mu1(k_ind+1,bk_ind+1,eps_j) = mu1(k_ind+1,bk_ind+1,eps_j) + (1-k_weight)*(1-bk_weight)*pie(eps_i,eps_j)*mu(k_i,bk_i,eps_i)*(1-exit);

                        % entrants
                        mu1(k_ent_ind,bk_ent_ind,eps_j) = ( mu1(k_ent_ind,bk_ent_ind,eps_j) + k_ent_weight*bk_ent_weight*pie(eps_i,eps_j)*mu(k_i,bk_i,eps_i)*exit );
                        mu1(k_ent_ind+1,bk_ent_ind,eps_j) = ( mu1(k_ent_ind+1,bk_ent_ind,eps_j) + (1-k_ent_weight)*bk_ent_weight*pie(eps_i,eps_j)*mu(k_i,bk_i,eps_i)*exit );
                        mu1(k_ent_ind,bk_ent_ind+1,eps_j) = ( mu1(k_ent_ind,bk_ent_ind+1,eps_j) + k_ent_weight*(1-bk_ent_weight)*pie(eps_i,eps_j)*mu(k_i,bk_i,eps_i)*exit );
                        mu1(k_ent_ind+1,bk_ent_ind+1,eps_j) = ( mu1(k_ent_ind+1,bk_ent_ind+1,eps_j) + (1-k_ent_weight)*(1-bk_ent_weight)*pie(eps_i,eps_j)*mu(k_i,bk_i,eps_i)*exit );

                    end
                end
            end
        end
    end



    dist_error = max(max(max(abs(mu1-mu))));
    mu = mu1;
end


% [Kmat, BKmat] = meshgrid(kgrid, bkgrid);
% medmu= v(1:knum,1:bknum,3);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

agg_K = 0;
agg_Y = 0;
agg_I = 0;
agg_N = 0;
agg_B = 0;
borrowing = zeros(knum, bknum, enum);
agg_Dive = borrowing;
total_bail = 0;
agg_Cash = 0;
agg_Dive = 0;


% labor
for k_i = 1:knum
    for bk_i = 1:bknum
        for eps_i = 1:enum
            agg_N = agg_N + n(k_i, bk_i, eps_i)*mu(k_i, bk_i, eps_i);
        end
    end
end


% cash
for k_i = 1:knum
    for bk_i = 1:bknum
        for eps_i = 1:enum
            agg_Cash = agg_Cash + cash(k_i, bk_i, eps_i)*mu(k_i, bk_i, eps_i);
        end
    end
end

% dive
for k_i = 1:knum
    for bk_i = 1:bknum
        for eps_i = 1:enum
            agg_Dive(k_i, bk_i, eps_i) = dive(k_i, bk_i, eps_i)*mu(k_i, bk_i, eps_i);
        end
    end
end



% bailout accounts

for k_i = 1:knum
    for bk_i = 1:bknum
        for eps_i = 1:enum
            total_bail = total_bail + bailout(k_i, bk_i, eps_i)*mu(k_i, bk_i, eps_i);
        end
    end
end
Gov = wss*agg_N*tax;
Gov_E = Gov - total_bail;

% output
for k_i = 1:knum
    for bk_i = 1:bknum
        for eps_i = 1:enum
            agg_Y = agg_Y + y_val(k_i, bk_i, eps_i)*mu(k_i, bk_i, eps_i);
        end
    end
end

% capital
for k_i = 1:knum
    agg_K = agg_K + kgrid(k_i)*sum(sum(mu(k_i,:,:)));
end

% for debt = leverage * capital
for k_i = 1:knum
    for bk_i = 1:bknum
         if bkgrid(bk_i) > 0
            agg_B = agg_B + (kgrid(k_i)*bkgrid(bk_i)) *sum(mu(k_i,bk_i,:));
         end
    end
end

% % for leverage
% for bk_i = 1:bknum
%     if bkgrid(bk_i) > 0
%         agg_B = agg_B + bkgrid(bk_i) *sum(sum(mu(:,bk_i,:)));
%     end
% end



% consumption
agg_C = agg_Y - agg_K + (1-delta)*agg_K - Gov_E;

% investment
agg_I = agg_Y - agg_C - Gov_E;
