% Create the Jacobian functions for computations in the Kalman Filter.
% Do not try to save functions to file, instead run this script at startup.

% Hard code feature points.

t = sym('t');
N = sym('N_',[3 4]);
N_c = sym('N_c_');
N_b = sym('N_b_',[3 1]);
xi = sym('xi_',[12 1]);
u = sym('u_',[6 1]);

F = jacobian(cntr_model_d(xi,u,tau),xi);
Fwru = jacobian(cntr_model_d(xi,u,tau),u);
H = jacobian(measure_model_d(t,xi,fp,f,N,N_c,N_b),...
    xi);
HwrN = jacobian(measure_model_d(t,xi,fp,f,N,N_c,N_b),...
    [N(:,1); N(:,2); N(:,3); N(:,4); N_c; N_b]);

disp('Generating symbolic functions depending on fp...')
disp('fp 1');disp(fp(1));
disp('fp 2');disp(fp(2));
disp('fp 3');disp(fp(3));
disp('fp 4');disp(fp(4));

Ffun = matlabFunction(F);
% function with handle @(xi_1,xi_2,xi_3,xi_8,xi_9,xi_10,xi_11,xi_12)

Fwrufun = matlabFunction(Fwru);

Hfun = matlabFunction(H);
% function with handle @(N_1_1,N_1_2,N_1_3,N_1_4...
%                        N_2_1,N_2_2,N_2_3,N_2_4...
%                        N_3_1,N_3_2,N_3_3,N_3_4...
%                        N_c,t,...
%                        xi_1,xi_2,xi_3,xi_7,xi_8,xi_9)

HwrNfun = matlabFunction(HwrN);
% function with handle @(N_1_1,N_1_2,N_1_3,N_1_4...
%                        N_2_1,N_2_2,N_2_3,N_2_4...
%                        N_3_1,N_3_2,N_3_3,N_3_4...
%                        N_c,t,...
%                        xi_1,xi_2,xi_3,xi_7,xi_8,xi_9)