clear all
close all
clc

dx0 = 0;
dy0 = 1;
dtheta0 = 0;
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
epsilon = 0.5;
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


tiles = tiledlayout(3,1);

nexttile(tiles);
line1 = plot(tlist_nonlinear, Vlist_nonlinear(:, 1)-Vlist_nonlinear(1,1),'DisplayName','Nonlinear');
hold on
line2 = plot(tlist_linear, Vlist_linear(:,1)-Vlist_nonlinear(1,1),'--','DisplayName','Linear');
xlabel('Time (s)')
ylabel('Position')
title('X Displacement')

nexttile(tiles);
plot(tlist_nonlinear, Vlist_nonlinear(:, 2)-Vlist_nonlinear(1,2))
hold on
plot(tlist_linear, Vlist_linear(:,2)-Vlist_nonlinear(1,2),'--')
xlabel('Time (s)')
ylabel('Position')
title('Y Displacement')

nexttile(tiles);
plot(tlist_nonlinear, Vlist_nonlinear(:, 3)-Vlist_nonlinear(1,3))
hold on
plot(tlist_linear, Vlist_linear(:,3)-Vlist_nonlinear(1,3),'--')
xlabel('Time (s)')
ylabel('Angle (rad)')
title('Rotation Angle')

tile_legend = legend([line1,line2]);
tile_legend.Layout.Tile = 'South';













