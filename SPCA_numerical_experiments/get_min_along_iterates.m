function G_min = get_min_along_iterates(G)
% G(k) = min(G(1:k)) --> best (smallest) (sub)gradient norm of F(x)
    N = length(G) - 1;
    G_min = NaN(1,N+1);
    G_min(1) = G(1);
    for i = 2 : N+1
        G_min(i) = min(G_min(i-1), G(i));
    end
end