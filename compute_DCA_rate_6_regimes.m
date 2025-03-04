function [p_active, id] = compute_DCA_rate_6_regimes(L1, mu1, L2, mu2)
id = 0; 
if mu1 + mu2 <= 0 && (mu1 ~= 0 && mu2 ~= 0) % No convergence
    p_active = 0;
    return;
else
    if (mu1 == 0) && (mu2 == 0)
        p_cvx = 1/L1 + 1/L2;
        p_active = p_cvx; id = 1;
        return;
    end
end

if L2 == inf
    p1 = (1-mu1/L2)/(L1-mu1) + 1/L2 * ( 1 + (1/L2-1/L1) / (1/mu1 - 1/L1) );
else
    p1 = 1/L2*(L2-mu1)/(L1-mu1) + 1/L2 * ( 1 + (1/L2-1/L1) / (1/mu1 - 1/L1) );
end

if L1 == inf
    p2 = (1-mu2/L1)/(L2-mu2) + 1/L1 * ( 1 + (1/L1-1/L2) / (1/mu2 - 1/L2) );
else
    p2 = 1/L1*(L1-mu2)/(L2-mu2) + 1/L1 * ( 1 + (1/L1-1/L2) / (1/mu2 - 1/L2) );
end

p3 = 1/(L2+mu2) + 1/L1 * ( 1 + (1/L1) / (1/mu2 + 1/L2 + (1/mu1 - 1/L1) ) );
p4 = mu1/mu2 * (1/mu2 + 1/mu1);
p5 = (mu1 + L2) / L2^2;
p6 = (mu2 + L1) / L1^2;

if mu1 == 0 || mu2 == 0
    thr_mu1_geq_0 = 0; thr_mu2_geq_0 = 0;
else
    thr_mu1_geq_0 = (1/mu1+1/mu2+1/L2);
    thr_mu2_geq_0 = (1/mu1+1/mu2+1/L1);
end

p_active = [];
if (L1 >= L2) && (L2 >= mu1) && (mu1 >= 0) && ...
      ( (mu2 >=0 ) || ( (mu2 < 0 ) && ( thr_mu1_geq_0 <= (2+L2/mu2)/L1 ) ) )
    p_active = p1;
    id = 1;
end
if (L2 >= L1) && (L1 > mu2) && (mu2 >= 0) && ...
      ( (mu1 >=0 ) || ( (mu1 < 0 ) && ( thr_mu2_geq_0 <= (2+L1/mu1)/L2 ) ) )
    p_active = p2;
    id = 2;
end
if (mu1 > -mu2) && (mu2 < 0) && (L2>mu1) && (L1>mu2) && ...
           ( thr_mu1_geq_0 > (2+L2/mu2)/L1 ) && ( thr_mu1_geq_0 <= 0 )
    p_active = p3;
    id = 3;
end
if (mu1 > -mu2) && (mu2 < 0) && (L1 >= mu2) && ...
    ( ( thr_mu1_geq_0 >= 0 && L2 > 0 ) || L2 <= 0 )
    p_active = p4;
    id = 4;
end
if (mu1 >= L2) && (L2 >= 0) && (mu2 * thr_mu1_geq_0 >= 0)
    p_active = p5;
    id = 5;
end
if (mu2 >= L1)
    p_active = p6;
    id = 6;
end

if id == 0 % no regime found, no convergence
    p_active = 0;
end

end