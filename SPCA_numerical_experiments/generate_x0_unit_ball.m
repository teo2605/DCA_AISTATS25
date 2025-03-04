function x0_unit_ball = generate_x0_unit_ball(n,num_d)
    % function to generate sample points within the Euclidean unit ball
    % n: dimension of x0
    % num_d: number of points
    rand_directions = randn(n,num_d);
    rand_directions = rand_directions ./ vecnorm(rand_directions,2,1);    
    rnd_radii = rand(1,num_d); % random radius
    x0_unit_ball = rand_directions .* rnd_radii;
end