function [eulzyx,S,rnd1,rnd2] = euler_angle_calc(xi,rnd1,rnd2)
% Inputs and outputs
% S               = S is a rotation matrix that matches 
%                   the condition S*xi(1:3) = [0;0;norm(xi(1:3))]
%                   Diagnostic in nature                                3x3
% eulzyx          = 3-2-1 euler angles to register as i.c. for Th       3x1
% rnd1, rnd2      = random inputs, in the case where center_att = 1, these
%                   are reused to compute S and thus eulzyx so new angles
%                   stay close

% Computes xyz euler angles based on r_a
    
    d = xi(1:3);
    eulzyx = [0; 0; 0];
    rnd_multipl = 1;
    
    switch nargin
        case 1
            while ~all(eulzyx) || ~isreal(eulzyx)
                try
                    A = zeros(3);
                    while rank(A) ~= 3 
                        rnd2 = normrnd(0,rnd_multipl,[3,1]);
                        rnd1 = normrnd(0,rnd_multipl,[3,1]);

                        A = [rnd1+rnd1.*d rnd2+rnd2.*d -d];
                        S = norm(d,2)*eye(3)/A;
                    end
                    
                    theta = asin(S(1,3));
                    phi = asin(-S(2,3)/cos(theta));
                    psi = asin(-S(1,2)/cos(theta));
                    
                    eulzyx = [phi; theta; psi];
                catch
                    eulzyx = [0; 0; 0];
                end
            end
            
        case 3
            A = [rnd1+rnd1.*d rnd2+rnd2.*d -d];
            S = norm(d,2)*eye(3)/A;
            
            theta = asin(S(1,3));
            phi = asin(-S(2,3)/cos(theta));
            psi = asin(-S(1,2)/cos(theta));

            eulzyx = [phi; theta; psi];
    end
    
end



    


