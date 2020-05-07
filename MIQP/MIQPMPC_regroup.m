%MIQP-MPC
% Works in positive z-direction. 

% S  <-  A,B_1,B_3  <-  A^d,C^0,B^d
% G  <-  E_1,E_2,E_3,E_4 <- minf,Minf,L,E,E^i,A^i  <-  C^1
% W  <-  E_5  <-  M^i  <-  constants
% T  <-  M  <-  constants
% Qbar  <-  tuning matrices
% Rbar  <-  tuning matrices

% Initialize properties
% ADJUSTABLE%%%%%%%%%%%%%%
s = 2;
tau = 5;
N = 20;

xi_len = 12+12+1+tau+12;
u_len = 6+1+s;

% Initialize data
accel_up_lim = 0.005; % units see paper
moment_up_lim = 0.05;
r_up_ground = 3;
r_up_ref = [11.75,11];
xi_up_bthresh = 10000;
r_up_h = r_up_ref;
r_up_l = r_up_ref-0.1;
m_inf = -1000000*ones(xi_len,1);
M_inf = 1000000*ones(xi_len,1);
m_inf1 = m_inf(1);
M_inf1 = M_inf(1);
m_inf3 = m_inf(1:3);
M_inf3 = M_inf(1:3);
M_small3 = 0.000001;

Q = diag([ones(12,1); 5*ones(12,1); 1; ones(tau,1); ones(12,1)]);
P = Q;
R = diag([0.5*ones(6,1); 50; ones(s,1)]);
R_1 = ones(s*xi_len);

xi_init = [5;-3;14;0;0;0;2.7405;-0.0001;-2.1169;0;0;0;...
    zeros(12,1);0;zeros(tau,1);5.5;-2.5;11.5;0;0;0;2.7405;-0.0001;-2.1169;0;0;0];
% END ADJUSTABLE%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize matrices 1st level
C_up_1 = zeros(tau,tau);
L = zeros(12+1+2*s,u_len);
E = zeros(12+1+2*s,xi_len);
M = zeros(5+2*s,1);
E_up = zeros(s,3,xi_len);
E_up_len = 3;
M_up = zeros(s,3,1);
A = zeros(xi_len);
B_1 = zeros(xi_len,u_len);
B_3 = zeros(xi_len,s*xi_len);
Atilde_up = zeros(s,xi_len,xi_len);
E_3 = zeros(size(E,1)+s*xi_len*4+s*E_up_len*2+2*s,s*xi_len);
E_1 = zeros(size(E_3,1),u_len);
E_4 = zeros(size(E_3,1),xi_len);
E_5 = zeros(size(E_3,1),1);

% C_up_1
for k = 1:tau-1
    C_up_1(k,k+1) = 1;
end

E(13:end,3) = -1;
E(14:2:end-1,3) = 1;

%L
for k = 1:3
    L(k,k) = 1;
    L(3+k,k) = -1;
    L(6+k,3+k) = 1;
    L(9+k,3+k) = -1;
end
L(12+1+1:end,7) = -1;
for k = 1:s
L(12+1+(k-1)*2+1:12+1+(k-1)*2+2,7+k) = [M_inf1; M_inf1];
end

%M
M(1:6) = accel_up_lim;
M(7:12) = moment_up_lim;
M(13) = r_up_ground;
for k = 1:s
M(13+2*(k-1)+1) = M_inf1+r_up_ref(k);
M(13+2*(k-1)+2) = M_inf1-r_up_ref(k);
end

%E_up
for k = 1:s
    E_up(k,1,12+12+1) = -1;
    E_up(k,2,3) = 1;
    E_up(k,3,3) = -1;
end

%M_up
for k = 1:s
    M_up(k,1,1) = xi_up_bthresh;
    M_up(k,2,1) = -r_up_h(k);
    M_up(k,3,1) = r_up_l(k);
end

%A
A(1:12,1:12) = A_up_d;

A(13:24,1:12) = eye(12);
A(13:24,end-11:end) = -eye(12);

A(24+1,1:12) = C_up_0;
A(24+1,24+1) = 1;

A(24+1+tau,1:12) = -C_up_0;
A(24+1+tau,24+1) = -1;

A(end-12+1:end,end-12+1:end) = eye(12);

%B_1
B_1(1:12,1:6) = B_up_d;

%B_3
for k = 1:s
    B_3(1:xi_len,(k-1)*xi_len+1:k*xi_len) = eye(xi_len);
end

%Atilde_up
for k = 1:s
    Atilde_up(k,24+1,1:12) = -C_up_0;
    Atilde_up(k,24+1,26) = 1;

    Atilde_up(k,26:26+tau-1,26:26+tau-1) = C_up_1;
    Atilde_up(k,24+1+tau,1:12) = C_up_0;
    Atilde_up(k,24+1+tau,24+1) = 1;
end

dsh1 = size(L,1);
dsh2 = dsh1+4*xi_len*s;

%E_3
for k = 1:s
E_3(dsh1+xi_len*(k-1)*4+1:dsh1+k*xi_len*4,xi_len*(k-1)+1:xi_len*k) = ...
    [eye(xi_len);-eye(xi_len);eye(xi_len);-eye(xi_len)];
end

%E_1
E_1(1:dsh1,:) = -L;
for k = 1:s
    E_1(dsh1+4*xi_len*(k-1)+1:dsh1+4*xi_len*k,7+k) = ...
        [M_inf-m_inf;M_inf-m_inf;m_inf-M_inf;m_inf-M_inf];
    E_1(dsh2+2*E_up_len*(k-1)+1:dsh2+2*E_up_len*k,7+k) = [-M_inf3;-m_inf3];
end
E_1(end-2*s+1:end-s,end-s+1:end) = eye(s);
E_1(end-s+1:end,end-s+1:end) = -eye(s);

%E_4
E_4(1:dsh1,:) = -E;
for k = 1:s
    E_4(dsh1+4*xi_len*(k-1)+1:dsh1+4*xi_len*k,:) = ...
        [zeros(xi_len);zeros(xi_len);squeeze(Atilde_up(k,:,:));squeeze(-Atilde_up(k,:,:))];
    E_4(dsh2+2*E_up_len*(k-1)+1:dsh2+2*E_up_len*k,:) = [squeeze(E_up(k,:,:));-squeeze(E_up(k,:,:))];
end

%E_5
E_5(1:dsh1) = M;
for k = 1:s
    E_5(dsh1+xi_len*(k-1)*4+1:dsh1+xi_len*k*4) = ...
        [zeros(2*xi_len,1);-(m_inf-M_inf);-(m_inf-M_inf)];
    E_5(dsh2+2*E_up_len*(k-1)+1:dsh2+2*E_up_len*k) = ...
        [M_inf3-squeeze(M_up(k,:,:))';M_inf3-squeeze(M_up(k,:,:))'];
end
E_5(end-2*s+1:end-s) = ones(s,1);
E_5(end-s+1:end) = ones(s,1);

%%%%%%%%%%%%%%%%%%%%%%%
% Initialize matrices second level

S = zeros(xi_len*N,(u_len+s*xi_len)*N);
F = zeros(xi_len*N,xi_len);
for k = 1:N 
    F(1:xi_len,:) = A^k;
     for kk = 1:k
         % First Third
         S(xi_len*(k-1)+1:xi_len*k,u_len*(kk-1)+1:u_len*kk) = A^(k-1-(kk-1))*B_1;
         S(xi_len*(k-1)+1:xi_len*k,...
             N*u_len+s*xi_len*(kk-1)+1:N*u_len+s*xi_len*kk) = ...
                A^(k-1-(kk-1))*B_3;
     end
end

G = zeros(size(E_3,1)*N,(size(E_1,2)+size(E_3,2))*N);
W = zeros(size(E_5,1)*N,size(E_5,2));
T = zeros(size(E_4,1)*N,size(F,2));

Gscndterm = zeros(size(G,1),N*size(E_4,2));
for k = 1:N 
     G((k-1)*size(E_1,1)+1:k*size(E_1,1),(k-1)*size(E_1,2)+1:k*size(E_1,2)) = -E_1;
     G((k-1)*size(E_3,1)+1:k*size(E_3,1),...
         size(E_1,2)*N+(k-1)*size(E_3,2)+1:size(E_1,2)*N+k*size(E_3,2)) =...
            E_3;
     Gscndterm((k-1)*size(E_4,1)+1:k*size(E_4,1),(k-1)*size(E_4,2)+1:k*size(E_4,2)) = E_4;
     
     W((k-1)*size(E_5,1)+1:k*size(E_5,1),:) = E_5;
               
end

G = G-Gscndterm*S;
T = Gscndterm*F;

%%%%%%%%%%%%%%%%%%%%%%%%

R_R_1 = zeros(size(R,1)+size(R_1,1));
R_R_1(1:size(R,1),1:size(R,2)) = R;
R_R_1(size(R,1)+1:end,size(R,2)+1:end) = R_1;

Qbar = zeros(size(Q,1)*N);
Rbar = zeros(size(R_R_1,1)*N);
for k = 1:N
    Qbar((k-1)*size(Q,1)+1:k*size(Q,1),(k-1)*size(Q,2)+1:k*size(Q,2)) = Q;
    Rbar((k-1)*size(R_R_1,1)+1:k*size(R_R_1,1),(k-1)*size(R_R_1,2)+1:k*size(R_R_1,2)) = R_R_1;
end


H = S'*Qbar*S+Rbar;
q = S'*Qbar*F*xi_init;


%%%%%%%%%%%%%%%%%%%%%%%
% Solve the QP
bigU = quadprog(2*H,2*q',G,W+T*xi_init);



%% Plotting
plot(bigU(7:(24+6+1+2):end/2),'k')
grid on
xlabel('Time Point During Descent');
ylabel('epsion')
title('History of Epsilon Over one Horizon')










                
    
































