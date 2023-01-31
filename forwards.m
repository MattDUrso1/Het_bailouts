% forwards
% needs     T, mu, yval, n, kf, bkf to work

function [agg_k, agg_b, agg_n, agg_y, Gov_Et, w_vec, q_vec, dist] = forwards(T, mu, y_val, n, kft, bkft, w_vec, q_vec, bailout_t);

[alpha, nu, delta, beta, qss, exit, psi, zeta, tax, opp_cost, zss, knum, kgrid, k_enter, bknum, bkgrid, bk_enter, enum, epsgrid, pie] = params();


%Initialize with steady state distribution
dist = mu;


agg_k = zeros(1,T+1); agg_b = agg_k; 
agg_y = zeros(1,T); agg_n = agg_y; Govt = agg_y; Gov_Et = agg_y; total_bail_t = agg_y;

for i_k = 1:knum
    agg_k(1) = agg_k(1) + kgrid(i_k)*sum(sum(dist(i_k,:,:)));
end
for i_b = 1:bknum
    agg_b(1) = agg_b(1) + bkgrid(i_b)*sum(sum(dist(:,i_b,:)));
end


for t = 1:T

    dist1 = zeros(knum,bknum,enum);

    for i_k = 1:knum
        for i_b = 1:bknum
            for i_eps = 1:enum
            	if dist(i_k,i_b,i_eps) > 0
            		%Note, we're using date-specific decision rules now
            		[j_k,vw_k] = linear(knum, kgrid, kft(i_k,i_b,i_eps,t));
            		[j_b,vw_b] = linear(bknum, bkgrid, bkft(i_k,i_b,i_eps,t));

                    % entrants
                    [k_ent_ind,k_ent_weight] = linear(knum, kgrid, k_enter);
                    [bk_ent_ind,bk_ent_weight] = linear(bknum, bkgrid, bk_enter);

            		for j_eps = 1:enum
            			dist1(j_k,j_b,j_eps) = dist1(j_k,j_b,j_eps) + vw_k*vw_b*pie(i_eps,j_eps)*dist(i_k,i_b,i_eps)*(1-exit);
            			dist1(j_k+1,j_b,j_eps) = dist1(j_k+1,j_b,j_eps) + (1-vw_k)*vw_b*pie(i_eps,j_eps)*dist(i_k,i_b,i_eps)*(1-exit);
            			dist1(j_k,j_b+1,j_eps) = dist1(j_k,j_b+1,j_eps) + vw_k*(1-vw_b)*pie(i_eps,j_eps)*dist(i_k,i_b,i_eps)*(1-exit);
            			dist1(j_k+1,j_b+1,j_eps) = dist1(j_k+1,j_b+1,j_eps) + (1-vw_k)*(1-vw_b)*pie(i_eps,j_eps)*dist(i_k,i_b,i_eps)*(1-exit);

                        % entrants
                        dist1(k_ent_ind,bk_ent_ind,j_eps) = ( dist1(k_ent_ind,bk_ent_ind,j_eps) + k_ent_weight*bk_ent_weight*pie(i_eps,j_eps)*dist(i_k,i_b,i_eps)*exit );
                        dist1(k_ent_ind+1,bk_ent_ind,j_eps) = ( dist1(k_ent_ind+1,bk_ent_ind,j_eps) + (1-k_ent_weight)*bk_ent_weight*pie(i_eps,j_eps)*dist(i_k,i_b,i_eps)*exit );
                        dist1(k_ent_ind,bk_ent_ind+1,j_eps) = ( dist1(k_ent_ind,bk_ent_ind+1,j_eps) + k_ent_weight*(1-bk_ent_weight)*pie(i_eps,j_eps)*dist(i_k,i_b,i_eps)*exit );
                        dist1(k_ent_ind+1,bk_ent_ind+1,j_eps) = ( dist1(k_ent_ind+1,bk_ent_ind+1,j_eps) + (1-k_ent_weight)*(1-bk_ent_weight)*pie(i_eps,j_eps)*dist(i_k,i_b,i_eps)*exit );
            		end
            	end
            end
        end
    end

    %Find aggregates
for k_i = 1:knum
    for bk_i = 1:bknum
        for eps_i = 1:enum
            agg_y(t) = agg_y(t) + y_val(k_i, bk_i, eps_i)*dist(k_i, bk_i, eps_i);
        end
    end
end

    for i_bk = 1:bknum
        for i_k = 1:knum
            for i_eps = 1:enum
            	agg_n(t) = agg_n(t) + n(i_k,i_bk,i_eps)*sum(dist(i_k,i_bk,i_eps));
            end
        end
    end

    for k_i = 1:knum
    for bk_i = 1:bknum
        for eps_i = 1:enum
            total_bail_t(t) = total_bail_t(t) + bailout_t(k_i, bk_i, eps_i, t)*mu(k_i, bk_i, eps_i);
        end
    end
    end
    Govt(t) = w_vec(t)*agg_n(t)*tax;
    Gov_Et(t) = Govt(t) - total_bail_t(t);

    %Future so dist1, not dist
    for i_k = 1:knum
    	agg_k(t+1) = agg_k(t+1) + kgrid(i_k)*sum(sum(dist1(i_k,:,:)));
    end

    for i_k = 1:knum
    for bk_i = 1:bknum
        if bkgrid(bk_i) > 0
    	agg_b(t+1) = agg_b(t+1) + (kgrid(i_k)*bkgrid(bk_i)) *sum(dist1(i_k,bk_i,:));
        end
    end
    end




    %Replace dist to start at next period
    dist = dist1;

%     c_guess(t) = agg_y(t) - agg_k(t+1) + (1-delta)*agg_k(t)
end
