function [ubar,xibar,Jval,feasible,output,U,Y,K,ref_adj] = control_alg(t,xi,xi_init,y,...
    Ffun,Fwrufun,Hfun,HwrNfun,P_NN,Xi_0,tau,betaa,A,H,V,center_att)

    % Selected Inputs and Outputs
    
    % tau           = sampling time                                     1x1
    % Fwrufun       = Function handle for Jacobian Evaluation of 
    %                 df/du. evaluates to B
    % betaa         = confidence percentage
    % Y_update      = determines if Y_tilde is updated at each iteration
    % A,H,V         = if not passed in evaluated at beginning of mission
    %                 time
    % center_att    = decision variable if attitude centering takes place
    
    % This book-keeps all actions necessary at each time in mission to compute u
    
    %- - - - - - - - - - - - - - - - - - - - - - - - - -                        
    num_noise = 16;
    num_foview = 16;
    num_horz = 20;
    num_x = 12;
    num_u = 6;
    %- - - - - - - - - - - - - - - - - - - - - - - - - -                        

    ref = evalin('base','ref');
    cnt = evalin('base','cnt');
    
    % Construct closed loop matricies based on current state
    B = Fwrufun();
   
    if A == 0
        A = Ffun(xi_init(1),xi_init(2),xi_init(3),...
                    xi_init(8),xi_init(9),xi_init(10),...
                    xi_init(11),xi_init(12));
    end   
    if H == 0
        H = Hfun(0,0,0,0,0,0,0,0,0,0,0,0,0,t,...
                   xi_init(1),xi_init(2),xi_init(3),...
                   xi_init(7),xi_init(8),xi_init(9));
        xi = xi_init;
    end   
    if V == 0
        V = HwrNfun(0,0,0,0,0,0,0,0,0,0,0,0,0,t,...
                   xi_init(1),xi_init(2),xi_init(3),...
                   xi_init(7),xi_init(8),xi_init(9));
    end
    
    G = 0.00001*ones(12,num_noise); %should represent model mismatch
    
    % Adjust the reference ------------------------------------------------
    % ---------------------------------------------------------------------
    
    if center_att
        % Adjust ref
        ref_adj = center_attitude(ref(:,cnt));
    else
        ref_adj = ref(:,cnt);
    end
        
    % Construct the initial-state-linearized QP ---------------------------
    % ---------------------------------------------------------------------
    
    % 1. Construct the tightened output constraints

        % Formulate additional terms in field of view constraint
        lin_area_xi = -H*xi+y;
        foview_subtract = zeros(num_foview,1);
        fov_cnt = 5;
        for i = 1:2:num_foview
            foview_subtract(i:i+1) = [lin_area_xi(fov_cnt);-lin_area_xi(fov_cnt)];
            fov_cnt = fov_cnt+1;
        end

        % Discretize
        sysNoise = ss(A,G,H,V,tau);
        sysControl = ss(A,B,zeros(1,12),0,tau);

        % Compute feedback gain
        K = place(sysControl.A,sysControl.B,...
            [0.99986;0.99987;0.99989;0.9999;0.99991;0.99992;...
            0.99993;0.99994;0.999906;0.99907;0.999908;0.9999]);

        % Evaluate Y_tilde_inf_beta
        [Y,U] = compute_constraints(1,foview_subtract,...
        sysControl.A-sysControl.B*K,sysNoise.B,sysNoise.C,sysNoise.D,K,...
        P_NN,Xi_0,num_horz,betaa);
        
    % 2. Formulate the QP and solve
        C = Y.A*[H;eye(num_x)];
        num_cnstr = size(Y.A,1);

        % Stack matrices only contains the initial A, C
        A0stack = zeros(num_x*num_horz,num_x);
        C0stack = zeros(num_cnstr*num_horz,num_x);
        xistack = zeros(num_x*num_horz,1);
        xirefstack = zeros(num_x*num_horz,1);
        urefstack = zeros(num_u*num_horz,1);
        
        for cnt = 1:num_horz
            A0stack(num_x*(cnt-1)+1:num_x*cnt,:) = A; 
            C0stack(num_cnstr*(cnt-1)+1:num_cnstr*cnt,:) = C;
            
            xirefstack(num_x*(cnt-1)+1:num_x*cnt) = xi;
            Bximatch = ref_adj-A0stack((cnt-1)*num_x+1:cnt*num_x,:)*ref_adj;
            urefstack((cnt-1)*num_u+1:cnt*num_u) = (B'*B)\B'*Bximatch;
        end
        
        for cnt = 1:num_horz
            xistack(num_x*(cnt-1)+1:num_x*cnt) = xi-ref_adj;
        end
        [ubar,xibar,Jval,feasible,output] = scvx_cnstr_slack_not_at_zero(...
            xistack,xirefstack,urefstack,A0stack,B,C0stack,Y,U,num_horz);  
        
end
