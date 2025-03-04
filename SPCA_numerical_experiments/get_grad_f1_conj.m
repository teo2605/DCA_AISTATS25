function df1_conj = get_grad_f1_conj(y,kappa,eta,lambda)
    if nargin < 4; lambda = 0; end
    N = length(y);
    if all(y == 0) || all(abs(y) == kappa)
    % in the case of second condition, the subdifferential is not a
    % singleton; it is {0} \cup {x \in B(0,1): sgn(xi) = sgn(yi), for all i}
        df1_conj = zeros(size(y)); 
    else
        get_pos_part = @(v) max(v,0); 
        df1_conj = ( sign(y) .* get_pos_part( abs(y) - kappa*ones(N,1) ) ) / ...
                   ( max( (eta-lambda), norm( get_pos_part( abs(y) - kappa*ones(N,1) ) ) ) );
    end
end