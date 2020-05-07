function y = measure_model_d(t,xi,fp,f,N_,N_c,N_b)
    % Inputs and outputs
    % xi            = [r v th Omega_b]'                                12x1
    % fp            = feature points [fp_1 ... fp_4]                    1x4 
    %                                              struct with fields x,y,z
    % f             = focal length of nav camera                        1x1
    % N                                                                 3x4
    % N_c                                                               1x1
    % N_b                                                               3x1
    % Y             = Measured outputs                                 15x1 
       
    % Constants
    omega=2*pi/(5.27*3600);        %rotation rate (rad/sec)
    
    if isa(xi,'sym')
        D_a = sym('D_a_',[3 4]);
        D_b = sym('D_b_',[3 4]);
        D_b_hat = sym('D_b_hat_',[3 4]);
        d_b_hat = sym('d_b_hat_',[4 1]);
        c = sym('c_',[2 4]);
    else
        D_a = zeros(3,4);
        D_b = zeros(3,4);
        D_b_hat = zeros(3,4);
        d_b_hat = zeros(4,1);
        c = zeros(2,4);
    end
    
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
        D_b_hat(:,k) = D_b(:,k)+(1-dot(-D_b(:,k)/norm(D_b(:,k),2),e_3))*N_(:,k);
        d_b_hat(k) = norm(D_b_hat(:,k),2);
        c(:,k) = ((f/([0 0 1]*D_b_hat(:,k)))+N_c)*[1 0 0;0 1 0]*D_b_hat(:,k);
    end
    
    Omega_b_hat = xi(10:12)+N_b;
    
    y = [d_b_hat;c(:,1);c(:,2);c(:,3);c(:,4);Omega_b_hat];
    

end


