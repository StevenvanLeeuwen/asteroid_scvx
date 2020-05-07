function [xi_t,F] = cntr_model_d(xi_tm1,u_tm1,tau)
    % Can handle either symbolic or numeric inputs
    % See notes for description of model
    
    % Constants
    mu=4.463e-04;                  %Gravitational constant for Eros (km^3/s^2)
    omega=2*pi/(5.27*3600);        %rotation rate (rad/sec)
    J = [500 0 0;0 500 0;0 0 500];       %interial matrix
    
%% Gravity Model Parameters

    x = xi_tm1(1);
    y = xi_tm1(2);
    z = xi_tm1(3);
    r=norm(xi_tm1(1:3));
    
    % ---------------------- Low Fidelity Gravity Model-----------------------
    % -------------------------------------------------------------------------
    
    a = 20; %km
    b = 5; %km
    c = 5; %km
    R = 20; %km
    C=zeros(5,5);
    S=zeros(5,5);
    

	% compute the intermediate parameters
    C(1,1) = 1;
	C(3,1) = 1/(5*R^2) * (c^2 - (a^2+b^2))/2;
	C(3,3) = (a^2 - b^2)/(20*R^2);
	C(5,1) = 15/7 * (C(3,1)^2 + 2*C(3,3)^2);
	C(5,3) = 5/7 * (C(3,1)*C(3,3));
	C(5,5) = 5/28*(C(3,3)^2);

% Montenbruck Recursions for Model

    %Follows the algorithm as laid out in O. Montenbruck and E. Gill: 
    % ‘‘Satellite Orbits: Models, Methods and Applications’’, 2011.

    if isa(xi_tm1,'sym')
        V = sym('V_',[6 6]);
        W = sym('W_',[6 6]);
        xdd = sym('xdd_',[6 6]);
        ydd = sym('xdd_',[6 6]);
        zdd = sym('xdd_',[6 6]);
        V(:,:) = 0;
        W(:,:) = 0;
        xdd(:,:) = 0;
        ydd(:,:) = 0;
        zdd(:,:) = 0;
    else
        V=zeros(6);
        W=zeros(6);
        xdd=zeros(5);
        ydd=zeros(5);
        zdd=zeros(5);
    end

    for j=1:6
        m = j-1;
        for i=j:6
            n = i-1;
            if j==i
                if i==1
                    V(i,j) = R/r;
                else
                    V(i,j)=(2*m-1)*((x*R/r^2)*V(i-1,j-1)-(y*R/r^2)*W(i-1,j-1));
                    W(i,j)=(2*m-1)*((x*R/r^2)*W(i-1,j-1)+(y*R/r^2)*V(i-1,j-1));
                end
            else
                if i==2
                    V(i,j)=((2*n-1)/(n-m))*(z*R/r^2)*V(i-1,j);
                    W(i,j)=((2*n-1)/(n-m))*(z*R/r^2)*W(i-1,j);
                else
                    V(i,j)=((2*n-1)/(n-m))*(z*R/r^2)*V(i-1,j)-((n+m-1)/(n-m))*(R^2/r^2)*V(i-2,j);
                    W(i,j)=((2*n-1)/(n-m))*(z*R/r^2)*W(i-1,j)-((n+m-1)/(n-m))*(R^2/r^2)*W(i-2,j);
                end
            end
        end
    end

% Gravity Partial Accelerations

    for i=1:5
        n=i-1;
        for j=1:i
            m=j-1;
            if m==0
                xdd(i,j)=mu/R^2*(-C(i,j)*V(i+1,2));
                ydd(i,j)=mu/R^2*(-C(i,j)*W(i+1,2));
            else
                xdd(i,j)=((mu/R^2)*0.5)*(-C(i,j)*V(i+1,j+1)-S(i,j)*W(i+1,j+1)+...
                    (factorial(n-m+2)/factorial(n-m))*(C(i,j)*V(i+1,j-1)+S(i,j)*W(i+1,j-1)));
                ydd(i,j)=((mu/R^2*0.5))*(-C(i,j)*W(i+1,j+1)+S(i,j)*V(i+1,j+1)+...
                    (factorial(n-m+2)/factorial(n-m))*(-C(i,j)*W(i+1,j-1)+S(i,j)*V(i+1,j-1)));
            end
            zdd(i,j)=((mu/R^2)*(n-m+1))*(-C(i,j)*V(i+1,j)-S(i,j)*W(i+1,j));
        end
    end
    gx=sum(sum(xdd));
    gy=sum(sum(ydd));
    gz=sum(sum(zdd));


    F=[gx;gy;gz];
    % ---------------------- End Low Fidelity Gravity Model-------------------
    % -------------------------------------------------------------------------
%% Copied from the last section of plant_c.m except for g.
    
    if isa(xi_tm1,'sym')
        dxi_tm1=sym('dxi_tm1_',[12 1]);
    else
        dxi_tm1 = zeros(12,1);
    end
    
    n = omega;
    Omegacap = [0;0;n];
    e_i = [1;0;0];
    e_j = [0;1;0];
    e_k = [0;0;1];
    
    R1i = [cos(xi_tm1(9)) -sin(xi_tm1(9)) 0;...
           sin(xi_tm1(9)) cos(xi_tm1(9)) 0;...
           0 0 1];
    R21 = [cos(xi_tm1(8)) 0 sin(xi_tm1(8));...
           0 1 0;...
           -sin(xi_tm1(8)) 0 cos(xi_tm1(8))];  

    B = [R21*R1i*e_i R1i*e_j e_k];
    
    dxi_tm1(1:3) = xi_tm1(4:6);
    dxi_tm1(4:6) = u_tm1(1:3)+F-2*cross(Omegacap,xi_tm1(4:6))-...
        cross(Omegacap,cross(Omegacap,xi_tm1(1:3)));
    dxi_tm1(7:9) = B\xi_tm1(10:12);
    dxi_tm1(10:12) = J\(u_tm1(4:6)-cross(xi_tm1(10:12),J*xi_tm1(10:12)));
    
    xi_t = xi_tm1 + tau.*dxi_tm1;

end

