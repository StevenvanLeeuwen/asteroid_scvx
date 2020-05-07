function [xi_return] = init_euler_calc(xi)
% Inputs and outputs
% xi               = [r_a ~ ~ ~]'                                       12x1
% xi_return        = [r_a ~ ~ ~]'                                       12x1
% rnd1, rnd2      = random inputs, in the case where center_att = 1, these
%                   are reused to compute S and eulzyx so new angles
%                   stay close

% This function continues to compute sets of Euler Angles that satisfy the
% condition (xi-avg_fp)(1:3) expressed in the body frame is [0;0;norm((xi-avg_fp)(1:3))],
% until they also satisfy the condition that Rba*xi(1:3) = [0;0;norm((xi-avg_fp)(1:3))]
% (with Rba = f(phi,theta,psi) ) with some tolerance

    fp = evalin('base','fp');
    avg_fp = mean([[fp(1).x;fp(2).x;fp(3).x;fp(4).x],...
         [fp(1).y;fp(2).y;fp(3).y;fp(4).y],...
         [fp(1).z;fp(2).z;fp(3).z;fp(4).z]])';
     
    while true

        [xi(7:9),~,rnd1,rnd2] = euler_angle_calc(xi(1:3)-avg_fp);
        assignin('base','rnd1',rnd1);
        assignin('base','rnd2',rnd2);
        
        omega = 0;
        t = 0;

        % Can change check condition as you like
        tol = norm(xi(1:3)-avg_fp,2)/40;

        % This is the convention of how our Euler Angles are defined.
        % Below must match that in init_focal_length.m
        Ria = [cos(omega*t) -sin(omega*t) 0;...
               sin(omega*t) cos(omega*t) 0;...
               0 0 1];
        R1i = [cos(xi(9)) -sin(xi(9)) 0;...
               sin(xi(9)) cos(xi(9)) 0;...
               0 0 1];
        R21 = [cos(xi(8)) 0 sin(xi(8));...
               0 1 0;...
               -sin(xi(8)) 0 cos(xi(8))];
        Rb2 = [1 0 0;...
               0 cos(xi(7)) -sin(xi(7));...
               0 sin(xi(7)) cos(xi(7))];

        Rba = Rb2*R21*R1i*Ria;

        % Want the transformation from eul to rotm to be accurate according to
        % tol
        check = Rba*(xi(1:3)-avg_fp);
        
        % Ideally, check should be [0;0;norm(xi(1:3))].
        if abs(check(1))<tol && abs(check(2))<tol && check(3) > 0
            xi_return = xi;
            break 
        end
        
    end

end
















