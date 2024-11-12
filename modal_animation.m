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
epsilon = 0.01;
tspan = [0 5];

my_rate_func = @(V_in) box_rate_func(tspan(1),V_in,box_params);
V0 = [dx0;dy0;dtheta0;vx0;vy0;vtheta0];
Veq = multi_newton_solver(my_rate_func,V0,true);
J_approx = approximate_jacobian(my_rate_func, Veq);

[Umode, omega_n] = eig(J_approx(4:6, 1:3));
Umode = Umode(:, 3);
omega_n = sqrt(-omega_n(1));
V0 = Veq + epsilon*[Umode;0;0;0];


tspan = [0 60];

%run the integration
my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
[tlist,Vlist] = ode45(my_rate_func,tspan,V0);

fname='modal_animation.avi';

% % create a videowriter, which will write frames to the animation file
% writerObj = VideoWriter(fname);
% 
% % must call open before writing any frames
% open(writerObj);
% 
% fig1 = figure(1);
box_corners = [-1 -1 1 1 -1; -1 1 1 -1 -1];
box_corners_world = compute_rbt(dx0, dy0, dtheta0, box_corners);
box = plot(box_corners_world(1,:),box_corners_world(2,:));
hold on;

num_zigs = 5;
w = .1;
num_springs = length(box_params.k_list);
for i = 1:num_springs
    spring_plots{i} = initialize_spring_plot(num_zigs,w);
end

axis equal; axis square;
axis([-3,3,-3,3]);
timestamp = text(-2.5, 2.5, '');
frametimes = zeros(length(tlist),1);

pause(2)
tic

for i = 1:length(tlist)
    t = tlist(i);
    V = Vlist(i,:)*10;
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
    % current_frame = getframe(fig1);
    % 
    % % write the frame to the video
    % writeVideo(writerObj,current_frame)
end

% close(writerObj);

% close all
% figure()
% plot(frametimes - tlist)
% xlabel('Timestep')
% ylabel('Frame Lag (s)')