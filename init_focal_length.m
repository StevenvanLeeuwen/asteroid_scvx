function f = init_focal_length(t,xi,fp)
    % Inputs and outputs
    % xi            = [r v rh omega_b]'                                12x1
    % fp            = feature points [fp_1 ... fp_4]                    1x4 
    %                                              struct with fields x,y,z
    % f             = focal length of nav camera                        1x1 
       
    
    % Choose among fp the distance to feature point vector that best projects
    % along the body z axis- then the body z distance to this feature point is
    % the init focal length.
    
    % Constants
    omega=2*pi/(5.27*3600);        %rotation rate (rad/sec)
    
    
    D_a = zeros(3,4);
    D_b = zeros(3,4);
    proj = zeros(4,1);
    
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
    
    e_3 = [0;0;1];

    for k = 1:4
        D_a(:,k) = [xi(1)-fp(k).x; xi(2)-fp(k).y; xi(3)-fp(k).z];
        D_b(:,k) = Rba*D_a(:,k);
        proj(k) = 1-dot(-D_b(:,k)/norm(D_b(:,k),2),e_3);
    end
    ind = min(proj)==proj;
    f = -D_b(3,ind);
    
end


