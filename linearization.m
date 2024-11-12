clear all
close all
clc

dx0 = 0;
dy0 = 1;
dtheta0 = -pi/6;
vx0 = 0;
vy0 = 0;
vtheta0 = 0;

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

%small number used to scale initial perturbation
epsilon = 0.005;
tspan = [0 30];

my_rate_func = @(V_in) box_rate_func(tspan(1),V_in,box_params);
V0 = [dx0;dy0;dtheta0;vx0;vy0;vtheta0];
Veq = multi_newton_solver(my_rate_func,V0,true);
V0 = Veq + epsilon*V0;

J_approx = approximate_jacobian(my_rate_func, V0);
my_linear_rate = @(t_in,V_in) J_approx*(V_in-Veq);

%run the integration of nonlinear system
% [tlist_nonlinear,Vlist_nonlinear] =...
% your_integrator(my_rate_func,tspan,V0,...);
my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
[tlist_nonlinear, Vlist_nonlinear] = ode45(my_rate_func,tspan,V0);

%run the integration of linear system
% [tlist_linear,Vlist_linear] =...
% your_integrator(my_linear_rate,tspan,V0,...);
[tlist_linear, Vlist_linear] = ode45(my_linear_rate,tspan,V0);

plot(tlist_nonlinear, Vlist_nonlinear(:, 1))
hold on
plot(tlist_linear, Vlist_linear(:, 1), '.--')