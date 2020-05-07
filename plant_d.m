function [xi_t,F] = plant_d(xi_tm1,u_tm1,tau)

% #codegen
% Will Dunham
% adopted from Eros_grav_rot_simu(x,y,z,vx,vy,vz,ux,uy,uz) in Dom's code
% only changed this section and past ending formulation of X(1:6)
    
    % See notes for description of model
    
    % Constants
    mu=4.463e-04;                  %Gravitational constant for Eros (km^3/s^2)
    omega=2*pi/(5.27*3600);        %rotation rate (rads/sec)
    J = [500 0 0;0 500 0;0 0 500];       %Inertial Matrix of satellite


%% Gravity Model Parameters

    x = xi_tm1(1);
    y = xi_tm1(2);
    z = xi_tm1(3);
    r=norm(xi_tm1(1:3));
    
    % ---------------------- High Fidelity Gravity Model-----------------------
    % -------------------------------------------------------------------------
    %Normalized Data from NEAR Mission
    %C(n,m)=factor(n,m)*C_norm(n,m)  
    %Same for S(n,m).
    %Data from Yeomans et al. 2001.  
    %Normalizing factors can be found at http://sebago.mit.edu/near/normal.html

    R=20;            %Reference radius for Eros (km)

    C=zeros(5,5);
    S=zeros(5,5);
    Pi=zeros(5,5);
    
    for i=1:5
        n=i-1;
        for j=1:i
            m=j-1;
            if m==0
                Pi(i,j)=sqrt(factorial(n+m)/(factorial(n-m)*(2*n+1)));
            else
                Pi(i,j)=sqrt(factorial(n+m)/(2*factorial(n-m)*(2*n+1)));
            end
        end
    end
    
    %Note that the degree i and order j would correspond to the (i+1,j+1) index
    %for C, S, and the Montenbruck Reccurssion terms V and W
    
    C(1,1)=1*1;                     S(1,1)=1*0;
    C(2,1)=1/Pi(2,1)*0;             S(2,1)=1/Pi(2,1)*0;
    C(2,2)=1/Pi(2,2)*0;             S(2,2)=1/Pi(2,2)*0;
    C(3,1)=1/Pi(3,1)*-0.052478;     S(3,1)=1/Pi(3,1)*0;
    C(3,2)=1/Pi(3,2)*0;             S(3,2)=1/Pi(3,2)*0;
    C(3,3)=1/Pi(3,3)*0.082533;       S(3,3)=1/Pi(3,3)*-0.027739;
    C(4,1)=1/Pi(4,1)*-0.001159;     S(4,1)=1/Pi(4,1)*0;
    C(4,2)=1/Pi(4,2)*0.004232;      S(4,2)=1/Pi(4,2)*0.003348;
    C(4,3)=1/Pi(4,3)*0.001834;      S(4,3)=1/Pi(4,3)*-0.000689;
    C(4,4)=1/Pi(4,4)*-0.010308;     S(4,4)=1/Pi(4,4)*-0.012218;
    C(5,1)=1/Pi(5,1)*0.012509;      S(5,1)=1/Pi(5,1)*0;
    C(5,2)=1/Pi(5,2)*-0.000105;     S(5,2)=1/Pi(5,2)*-0.00005;
    C(5,3)=1/Pi(5,3)*-0.017488;     S(5,3)=1/Pi(5,3)*0.004872;
    C(5,4)=1/Pi(5,4)*0.000056;      S(5,4)=1/Pi(5,4)*-0.000332;
    C(5,5)=1/Pi(5,5)*0.017534;      S(5,5)=1/Pi(5,5)*-0.008993;


% Montenbruck Recursions for Model

    %Follows the algorithm as laid out in O. Montenbruck and E. Gill: 
    % ‘‘Satellite Orbits: Models, Methods and Applications’’, 2011.

    V=zeros(6);
    W=zeros(6);

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
    xdd=zeros(5);
    ydd=zeros(5);
    zdd=zeros(5);

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

    F = [gx;gy;gz];
    % ---------------------- End High Fidelity Gravity Model-------------------
    % -------------------------------------------------------------------------
%% Finish rest of model, consistent with prediction except for gravity

    dxi_tm1 = zeros(12,1);
    
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
   
