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
epsilon = 1e-7;
tspan = [0 5];

my_rate_func = @(V_in) box_rate_func(tspan(1),V_in,box_params);
V0 = [dx0;dy0;dtheta0;vx0;vy0;vtheta0];

Veq = multi_newton_solver(my_rate_func,V0,true);
J_approx = approximate_jacobian(my_rate_func, Veq);

mode_index = 2;
[Umode, omega_n] = eig(J_approx(4:6, 1:3));
Umode = Umode(:, mode_index);
omega_n = sqrt(-omega_n(mode_index,mode_index));
V0 = Veq + epsilon*[Umode;0;0;0];

%run the integration of nonlinear system
% [tlist_nonlinear,Vlist_nonlinear] =...
% your_integrator(my_rate_func,tspan,V0,...);
my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
[tlist_nonlinear, Vlist_nonlinear] = ode45(my_rate_func,tspan,V0);

x_modal = Veq(1)+epsilon*Umode(1)*cos(omega_n*tlist_nonlinear);
y_modal = Veq(2)+epsilon*Umode(2)*cos(omega_n*tlist_nonlinear);
theta_modal = Veq(3)+epsilon*Umode(3)*cos(omega_n*tlist_nonlinear);

tiles = tiledlayout(3,1);

nexttile(tiles);
line1 = plot(tlist_nonlinear, Vlist_nonlinear(:, 1)-x_modal(1),'DisplayName','Nonlinear');
hold on
line2 = plot(tlist_nonlinear, x_modal-x_modal(1),'DisplayName','Modal');
xlabel('Time (s)')
ylabel('Position')
title('X Displacement')

nexttile(tiles);
plot(tlist_nonlinear, Vlist_nonlinear(:, 2)-y_modal(1))
hold on
plot(tlist_nonlinear, y_modal-y_modal(1))
xlabel('Time (s)')
ylabel('Position')
title('Y Displacement')

nexttile(tiles);
plot(tlist_nonlinear, Vlist_nonlinear(:, 3)-theta_modal(1))
hold on
plot(tlist_nonlinear, theta_modal-theta_modal(1))
xlabel('Time (s)')
ylabel('Angle (rad)')
title('Rotation Angle')

tile_legend = legend([line1,line2]);
tile_legend.Layout.Tile = 'South';