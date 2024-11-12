%Implementation of finite difference approximation
%for Jacobian of multidimensional function
%INPUTS:
%fun: the mathetmatical function we want to differentiate
%x: the input value of fun that we want to compute the derivative at
%OUTPUTS:
%J: approximation of Jacobian of fun at x
function J = approximate_jacobian(fun,x)
    size_x = length(x);
    size_y = length(fun(x));
    J = zeros([size_y,size_x]);
    h = 1e-6; % Kevin's rec
    for j = 1:size_x
        e_j = zeros([size_x,1]);
        e_j(j) = 1;
        J(:,j) = (fun(x + e_j*h) - fun(x - e_j*h)) ./ (2*h);
    end
end