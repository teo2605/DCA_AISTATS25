function [x] = linearoptimization_step(dir,ind,tag)
%
% [x] = linearoptimization_step(dir,ind)
%
% This routine performs a linear optimization step with objective function
% given by dir*x on the indicator function ind. That is, it evaluates
%   x=argmin_{ind(x)=0} [dir*x].
%
% Input: - direction dir, (gradient of the linear objective function)
%        - indicator function ind
%        - optional tag
%
% Output: x=argmin_{ind(x)=0} [dir*x].
%
g_imp=-dir;
x=Point('Point');
feas=Point('Function value');
ind.AddComponent(x,g_imp,feas);

if nargin > 2
    ind.AddComponent(x,g_imp,feas,tag);
else
    ind.AddComponent(x,g_imp,feas);
end
end

