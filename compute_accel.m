%Computes the linear and angular acceleration of the box
%given its current position and orientation
%INPUTS:
%x: current x position of the box
%y: current y position of the box
%theta: current orientation of the box
%box_params: a struct containing the parameters that describe the system
%Fields:
%box_params.m: mass of the box
%box_params.I: moment of inertia w/respect to centroid
%box_params.g: acceleration due to gravity
%box_params.k_list: list of spring stiffnesses
%box_params.l0_list: list of spring natural lengths
%box_params.P_world: 2 x n list of static mounting
% points for the spring (in the world frame)
%box_params.P_box: 2 x n list of mounting points
% for the spring (in the box frame)
%
%OUTPUTS
%ax: x acceleration of the box
%ay: y acceleration of the box
%atheta: angular acceleration of the box
function [ax,ay,atheta] = compute_accel(x,y,theta,box_params)
    P_world = box_params.P_world;
    P_box = box_params.P_box;
    k = box_params.k_list;
    l0 = box_params.l0_list;
    m = box_params.m;
    g = box_params.g;
    I = box_params.I;

    % Box mounting points transformed to world frame
    P_box_world = compute_rbt(x,y,theta,P_box);

    % Sum of forces
    F = compute_spring_force(k,l0,P_world,P_box_world);

    % XY acceleration
    F_ext = sum(F,2) + [0;-m*g];
    ax = F_ext(1) / m;
    ay = F_ext(2) / m;

    % Angular acceleration
    r = P_box_world - [x;y];
    r_pad = [r;zeros(1,size(r,2))];
    F_pad = [F;zeros(1,size(F,2))];
    T = norm(cross(r_pad,F_pad));
    atheta = T/I;
end











