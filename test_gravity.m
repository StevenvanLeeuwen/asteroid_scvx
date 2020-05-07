% test gravity model 
% Constants
    mu=4.463e-04;                  %Gravitational constant for Eros (km^3/s^2)
    
    x = -1.05;
	y = 0.85;
	z = 5.62;
	xdot = 0;
	ydot = 0;
    r=norm([x y z]);
    
    %%
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

    Fhighfid =[gx;gy;gz]
    % ---------------------- End High Fidelity Gravity Model-------------------
    % -------------------------------------------------------------------------
    %%
    % ---------------------- Low Fidelity Gravity Model-----------------------
    % ------------------------------------------------------------------------
    a = 20; %km
    b = 5; %km
    c = 5; %km
    R = 20; %km
    C=zeros(5,5);
    S=zeros(5,5);
    Pi=zeros(5,5);
    

	% compute the intermediate parameters
	r = sqrt(x^2 + y^2 + z^2);
    C(1,1) = 1;
	C(3,1) = 1/(5*R^2) * (c^2 - (a^2+b^2))/2;
	C(3,3) = (a^2 - b^2)/(20*R^2);
	C(5,1) = 15/7 * (C(3,1)^2 + 2*C(3,3)^2);
	C(5,3) = 5/7 * (C(3,1)*C(3,3));
	C(5,5) = 5/28*(C(3,3)^2);

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

    Flowfid=[gx;gy;gz]

    % ---------------------- End Low Fidelity Gravity Model-------------------
    % ------------------------------------------------------------------------

%% Low Fidelity Gravity Model (Controller Model) in NMPC Code

% gravity_force([x,y,z],mu,a,b,c,R);
% 
% function F3 = gravity_force(xi,mu,a,b,c,R)
% 
% 	x = xi(1);
% 	y = xi(2);
% 	z = xi(3);
% 
% 	% compute the intermediate parameters
% 	r = sqrt(x^2 + y^2 + z^2);
% 	g = c;
% 	C20 = 1/(5*R^2) * (c^2 - (a^2+b^2))/2;
% 	C22 = (a^2 - b^2)/(20*R^2);
% 	C40 = 15/7 * (C20^2 + 2*C22^2);
% 	C42 = 5/7 * (C20*C22);
% 	C44 = 5/28*(C22^2);
% 
%     F3 = zeros(3,1);
%     
% 	F3(1) = (mu*((C22*((30*R^4*x*y^2)/r^7 - (5*R*x*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/2 + (C44*((9*R*x*((7*R*x*((30*R^4*x*y^2)/r^7 - (5*R*x*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/r^2 + (7*R*y*((30*R^4*x^2*y)/r^7 + (5*R*y*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/r^2))/r^2 + (9*R*y*((7*R*x*((30*R^4*x^2*y)/r^7 + (5*R*y*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/r^2 - (7*R*y*((30*R^4*x*y^2)/r^7 - (5*R*x*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/r^2))/r^2))/2 + C44*((7*R^2*((30*R^4*x*y^2)/r^7 - (5*R*x*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/(2*r^2) - (63*R^2*z^2*((30*R^4*x*y^2)/r^7 - (5*R*x*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/(2*r^4)) + C20*((3*R^4*x)/(2*r^5) - (15*R^4*x*z^2)/(2*r^7)) - C22*((3*R^4*x)/(2*r^5) - (15*R^4*x*z^2)/(2*r^7)) - C40*((5*R^2*((3*R^4*x)/(2*r^5) - (15*R^4*x*z^2)/(2*r^7)))/(4*r^2) - (9*R*z*((4*R^5*x*z)/r^7 + (7*R*z*((3*R^4*x)/(2*r^5) - (15*R^4*x*z^2)/(2*r^7)))/(3*r^2)))/(4*r^2)) - (R^2*x)/r^3 - (3*(a^2 - b^2)*((5*R^2*((3*R^4*x)/(2*r^5) - (15*R^4*x*z^2)/(2*r^7)))/(4*r^2) - (9*R*z*((4*R^5*x*z)/r^7 + (7*R*z*((3*R^4*x)/(2*r^5) - (15*R^4*x*z^2)/(2*r^7)))/(3*r^2)))/(4*r^2))*(a^2/2 + b^2/2 - g^2))/(70*R^4) + ((a^2 - b^2)*((7*R^2*((30*R^4*x*y^2)/r^7 - (5*R*x*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/(2*r^2) - (63*R^2*z^2*((30*R^4*x*y^2)/r^7 - (5*R*x*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/(2*r^4))*(a^2/2 + b^2/2 - g^2))/(280*R^4)))/R^2;
% 
% 	F3(2) = (mu*((C44*((9*R*x*((7*R*x*((30*R^4*x*y^2)/r^7 - (5*R*x*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/r^2 + (7*R*y*((30*R^4*x^2*y)/r^7 + (5*R*y*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/r^2))/r^2 + (9*R*y*((7*R*x*((30*R^4*x^2*y)/r^7 + (5*R*y*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/r^2 - (7*R*y*((30*R^4*x*y^2)/r^7 - (5*R*x*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/r^2))/r^2))/2 - (C22*((30*R^4*x^2*y)/r^7 + (5*R*y*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/2 - C44*((7*R^2*((30*R^4*x*y^2)/r^7 - (5*R*x*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/(2*r^2) - (63*R^2*z^2*((30*R^4*x*y^2)/r^7 - (5*R*x*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/(2*r^4)) + C20*((3*R^4*y)/(2*r^5) - (15*R^4*y*z^2)/(2*r^7)) + C22*((3*R^4*y)/(2*r^5) - (15*R^4*y*z^2)/(2*r^7)) - (R^2*y)/r^3 - (5*C40*R^2*((3*R^4*y)/(2*r^5) - (15*R^4*y*z^2)/(2*r^7)))/(4*r^2) + (3*(a^2 - b^2)*((5*R^2*((3*R^4*x)/(2*r^5) - (15*R^4*x*z^2)/(2*r^7)))/(4*r^2) - (9*R*z*((4*R^5*x*z)/r^7 + (7*R*z*((3*R^4*x)/(2*r^5) - (15*R^4*x*z^2)/(2*r^7)))/(3*r^2)))/(4*r^2))*(a^2/2 + b^2/2 - g^2))/(70*R^4) + ((a^2 - b^2)*((7*R^2*((30*R^4*x*y^2)/r^7 - (5*R*x*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/(2*r^2) - (63*R^2*z^2*((30*R^4*x*y^2)/r^7 - (5*R*x*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/(2*r^4))*(a^2/2 + b^2/2 - g^2))/(280*R^4)))/R^2; 
% 
% 	F3(3) = (mu*(5*C40*((4*R^2*((2*R^4*z)/(3*r^5) + (5*R*z*(R^3/(2*r^3) - (3*R^3*z^2)/(2*r^5)))/(3*r^2)))/(5*r^2) + (9*R*z*((3*R^2*(R^3/(2*r^3) - (3*R^3*z^2)/(2*r^5)))/(4*r^2) - (7*R*z*((2*R^4*z)/(3*r^5) + (5*R*z*(R^3/(2*r^3) - (3*R^3*z^2)/(2*r^5)))/(3*r^2)))/(4*r^2)))/(5*r^2)) - 3*C20*((2*R^4*z)/(3*r^5) + (5*R*z*(R^3/(2*r^3) - (3*R^3*z^2)/(2*r^5)))/(3*r^2)) + (R^2*z)/r^3 + (3*((10*R^3*z*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^4 + (3*R*z*((5*R^2*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/(2*r^2) - (35*R^2*z^2*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/(2*r^4)))/r^2)*(a^2 - b^2)*(a^2/2 + b^2/2 - g^2))/(140*R^4) + (5*C22*R*z*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2 - (9*C44*R*z*((7*R*x*((30*R^4*x*y^2)/r^7 - (5*R*x*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/r^2 + (7*R*y*((30*R^4*x^2*y)/r^7 + (5*R*y*((3*R^3*x^2)/r^5 - (3*R^3*y^2)/r^5))/r^2))/r^2))/r^2))/R^2;
%     F3
% end

