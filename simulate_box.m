clear all
close all
clc

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

%load the system parameters into the rate function
%via an anonymous function
my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);

x0 = 0;
y0 = 1;
theta0 = 0;
vx0 = 0;
vy0 = 0;
omega0 = 0;

V0 = [x0;y0;theta0;vx0;vy0;omega0];

my_rate_func1 = @(V_in) box_rate_func(zeros(length(V_in)),V_in,box_params);

V0 = multi_newton_solver(my_rate_func1,V0,true);


tspan = [0 60];

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
timestamp = text(-4, 4, '');
frametimes = zeros(length(tlist),1);

    t = tlist(1);
    V = Vlist(1,:);
    box_corners_world = compute_rbt(V(1), V(2), V(3), box_corners);
    set(box,'XData',box_corners_world(1,:),'YData',box_corners_world(2,:))
    
    P_box_world = compute_rbt(V(1), V(2), V(3),box_params.P_box);
    for j = 1:num_springs
        update_spring_plot(spring_plots{j},P_box_world(:,j),box_params.P_world(:,j));
    end
    drawnow
pause(2)
tic

for i = 1:length(tlist)
    t = tlist(i);
    V = Vlist(i,:);
    box_corners_world = compute_rbt(V(1), V(2), V(3), box_corners);
    set(box,'XData',box_corners_world(1,:),'YData',box_corners_world(2,:))
    
    P_box_world = compute_rbt(V(1), V(2), V(3),box_params.P_box);
    for j = 1:num_springs
        update_spring_plot(spring_plots{j},P_box_world(:,j),box_params.P_world(:,j));
    end

    set(timestamp,'String',string(round(t,2)));

    frametimes(i) = double(toc);
    pause(t - frametimes(i));

    drawnow
end

close all
figure()
plot(frametimes - tlist)
xlabel('Timestep')
ylabel('Frame Lag (s)')

