%computes the force exerted by the spring at one of its ends
%INPUTS:
%PA: the position of the first end of the spring
%PB: the position of the second end of the spring
%k: the stiffness of the spring
%l0: the natural length of the spring
%OUTPUTS:
%F: the force exerted by the spring at end B
function F = compute_spring_force(k,l0,PA,PB)
    springs = length(k);
    l = zeros(1,springs);
    e_s = zeros(2,springs);
    F = zeros(2,springs);
    for i = 1:springs
        %current length of the spring
        l(i) = norm(PB(:,i)-PA(:,i));
        %unit vector pointing from PA to PB
        e_s(:,i) = (PB(:,i)-PA(:,i))/l(i);
        %Force exerted by spring at point B
        F(:,i) = -k(i)*(l(i)-l0(i))*e_s(:,i);
    end
end