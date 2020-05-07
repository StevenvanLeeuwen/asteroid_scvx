%% Initialize 
addpath('./eros_setup');
disp('Initializing...')

% load feature points
load fp_close.mat

% Constants
omega=2*pi/(5.27*3600);        %rotation rate (rads/sec)
    
% time
tfinal = 250; %s
tau = 1;
tspan = 0:tau:tfinal;
tdescent = 200;

% state
xi_init = [-1.57;1.32;6.75;0;0;0;0;0;0;0;0;0]; %km km/s rad rad/s
t = 0;
xi = xi_init;
ref_adj = xi_init;

% Calculate the initial attitude
[xi_init] = init_euler_calc(xi_init);

% Reference building
ref = zeros(12,length(tspan));
cnt = 1;
for t = tspan
    
    Ria = [cos(omega*t) -sin(omega*t) 0;...
           sin(omega*t) cos(omega*t) 0;...
           0 0 1];
    
    ref(:,cnt) = xi_init;
    th_ref = Ria*xi_init(7:9);
    if t <= tdescent
        ref(:,cnt) = [-1.07+((tdescent-t)/tdescent)*(xi_init(1)--1.07);...
            0.82+((tdescent-t)/tdescent)*(xi_init(2)-0.82);...
            5.8+((tdescent-t)/tdescent)*(xi_init(3)-5.8);...
            0;0;0;th_ref;0;0;0];
    else
        ref(:,cnt) = [-1.07;...
            0.82;...
            5.8;...
            0;0;0;th_ref;0;0;0];
    end
    cnt = cnt+1;
end

% Starting focal length is z of D_b(:,3) at time 0
f = init_focal_length(0,xi_init,fp);

% Get Jacobian function handles
gen_sym();
clearvars -except Ffun Fwrufun Hfun HwrNfun f fp xi_init tfinal omega ...
           ref tspan tau rnd1 rnd2

% Probability value within constraints
betaa = 0.95;

% Measurement Noise standard deviation and correlation
Sigma = [0.001*ones(12,1);0.001;0.001.*ones(3,1)];
Cor = eye(16);

%% Starting Covariance of State Error
Sig = ekf_init();

%% Simulate
disp('Simulating...')

xi = xi_init;
xi_nonvis = xi_init;
xi_record = zeros(length(tspan),12);
xi_nonvis_record = zeros(length(tspan),12);
u_record = zeros(length(tspan),6);
y_record = zeros(length(tspan),15);
y_fromref_record = zeros(length(tspan),15);
y_guess_record = zeros(length(tspan),15);
Sig_record = zeros(length(tspan),12,12);
ref_adj_record = zeros(12,length(tspan));

k = 1;
cnt = 1;
u = zeros(6,1);

ref_adj = xi_init;

for t = tspan
    
    if mod(k,20) == 0
        disp(['Iteration:',string(k),'at time:',string(tspan(k)),'sec']);
    end
    
    % 1. The spacecraft moves by influence of the asteroid gravity and
    %    thrusters
        [xi_nonvis,g] = plant_d(xi_nonvis,u,tau);

    % 2. Measurements are taken
        % 2.a. Generate measurement noise (Note P_NN doesn't change)
        [N_,N_c,N_b,P_NN] = measurement_noise_gen(Sigma,Cor);
        % 2.b Record measurements
        y = measure_model_d(t,xi_nonvis,fp,f,N_,N_c,N_b);
        y_fromref = measure_model_d(t,ref_adj,fp,f,0*N_,0*N_c,0*N_b);

    % 3. Obtain State Estimate
        [xi,Sig,y_guess,A,H,V] = ekf_naddnoise_d(...
              t,y,u,xi,Sig,P_NN,Ffun,Hfun,HwrNfun,fp,f,tau);
               
    % 4. Compute u by MPC law
        
        Y_update = 1;
        center_att = 1;

        [ubar,xibar,Jval,feasible,output,U,Y,K,ref_adj] = control_alg(t,xi,xi_init,y_guess,...
            Ffun,Fwrufun,Hfun,HwrNfun,P_NN,Sig,tau,betaa,A,H,V,center_att);
        u = ubar(1:6);

        
        if ~feasible
            disp('    ERROR: LQMPC Not feasible')
            tspan = tspan(1:k-1);
            ref = ref(:,1:k-1);
            xi_record = xi_record(1:k-1,:);
            xi_nonvis_record = xi_nonvis_record(1:k-1,:);
            y_record = y_record(1:k-1,:);
            y_fromref_record = y_fromref_record(1:k-1,:);
            y_guess_record = y_guess_record(1:k-1,:);
            u_record = u_record(1:k-1,:);
            Sig_record = Sig_record(1:k-1,:,:);
            ref_adj_record = ref_adj_record(:,1:k-1);
            return
        end
        
    % Record things
    xi_record(k,:) = xi';
    xi_nonvis_record(k,:) = xi_nonvis';
    y_record(k,:) = y';
    y_fromref_record(k,:) = y_fromref';
    y_guess_record(k,:) = y_guess';
    u_record(k,:) = u';
    Sig_record(k,:,:) = Sig;
    ref_adj_record(:,k) = ref_adj;
    k = k+1;
    cnt = cnt+1;
end

%% Plot
plot_all
%%
















