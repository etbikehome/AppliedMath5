function find_equilibrium()


x0 = 0;
y0 = 1;
theta0 = -pi/6;
vx0 = 0;
vy0 = 0;
omega0 = 0;

x_guess = [x0;y0;theta0;vx0;vy0;omega0];

%define system parameters
box_params = struct();
box_params.m = 1;
box_params.g = 1;
box_params.I = 4;
box_params.k_list = [1 1 1 1] * 4;
box_params.l0_list = [1 1 1 1];
box_params.P_world = [-2 2 -2 2;
                    -2 -2 2 2];
box_params.P_box = [-1 1 -1 1;
                    -1 -1 1 1];

my_rate_func = @(V_in) box_rate_func(zeros(length(V_in)),V_in,box_params);

multi_newton_solver(my_rate_func,x_guess,true)


end