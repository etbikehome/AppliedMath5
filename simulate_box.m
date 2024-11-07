clear all
close all
clc

%define system parameters
box_params = struct();
box_params.m = 1;
box_params.g = 1;
box_params.I = 1;
box_params.k_list = [2 2 2 2];
box_params.l0_list = [1 1 1 1];
box_params.P_world = [-2 2 -2 2;
                    -2 -2 2 2];
box_params.P_box = [-1 1 -1 1;
                    -1 -1 1 1];

%load the system parameters into the rate function
%via an anonymous function
my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);

x0 = 0;
y0 = 1;
theta0 = -pi/6;
vx0 = 0;
vy0 = 0;
omega0 = 0;

V0 = [x0;y0;theta0;vx0;vy0;omega0];
tspan = [0 50];

%run the integration
[tlist,Vlist] = ode45(my_rate_func,tspan,V0);

box_corners = [-1 -1 1 1 -1; -1 1 1 -1 -1];
box_corners_world = compute_rbt(x0, y0, theta0, box_corners);
box = plot(box_corners_world(1,:),box_corners_world(2,:));
hold on;

num_zigs = 5;
w = .1;
num_springs = length(box_params.k_list);
for i = 1:num_springs
    spring_plots{i} = initialize_spring_plot(num_zigs,w);
end

axis equal; axis square;
axis([-5,5,-5,5]);

for i = 1:length(tlist)
    tic
    t = tlist(i);
    V = Vlist(i,:);
    box_corners_world = compute_rbt(V(1), V(2), V(3), box_corners);
    set(box,'XData',box_corners_world(1,:),'YData',box_corners_world(2,:))
    
    P_box_world = compute_rbt(V(1), V(2), V(3),box_params.P_box);
    for i = 1:num_springs
        update_spring_plot(spring_plots{i},P_box_world(:,i),box_params.P_world(:,i));
    end

    drawnow
    pause(1/15 - toc);
end


