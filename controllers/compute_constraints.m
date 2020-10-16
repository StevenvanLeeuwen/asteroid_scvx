function [Y_tilde_inf_beta,U] = compute_constraints(Y_update,foview_subtract,...
                                A,G,H,V,K,...
                                P_NN,Xi_0,k_prop,betaa)

    %- - - - - - - - - - - - - - - - - - - - - - - - - -                        
    % Y 
    %   1+4*4=17 constraints on 15+12=27 outputs
    num_cnstr = 17;
    num_output = 27;
    num_measure = 15;
    num_foview = 16;
    % U
    %   12 constraints on 6 channels
    num_u = 6;
    num_u_cnstr = 12;
    %- - - - - - - - - - - - - - - - - - - - - - - - - -                        
    
    %%%%%%%%%% Constraint parameters to modify %%%%%%%%%%
    z_cnstr = 5.61;
    radius_cnstr = 0.3;
    u_trans_cnstr = 0.002;
    u_rot_cnstr = 0.5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    assignin('base','z_cnstr',z_cnstr);
    assignin('base','radius_cnstr',radius_cnstr);
    assignin('base','u_trans_cnstr',u_trans_cnstr);
    assignin('base','u_rot_cnstr',u_rot_cnstr);

    if Y_update == 1
        
        % Formulate Y -----------------------------------------------------
        % -----------------------------------------------------------------
        
        % initialize
        S = zeros(num_cnstr,num_output);
        s_ = zeros(num_cnstr,1); % deterministic
        s_k = zeros(num_cnstr,1);    

        % field of view constraint
        s_(1:num_foview) = radius_cnstr*ones(num_foview,1)-foview_subtract;
        S(1:4,5:6) = [1 0; -1 0; 0 1; 0 -1];
        S(5:8,7:8) = [1 0; -1 0; 0 1; 0 -1];
        S(9:12,9:10) = [1 0; -1 0; 0 1; 0 -1];
        S(13:16,11:12) = [1 0; -1 0; 0 1; 0 -1];
        
        % z constraint
        S(num_foview+1,num_measure+3) = -1;
        s_(num_foview+1) = -z_cnstr;
        Y = Polyhedron('A',S,'B',s_);
        
        % Forumulate U-----------------------------------------------------
        %------------------------------------------------------------------
        
        % initialize
        u_s_ = zeros(num_u_cnstr,1);
        u_s_k = zeros(num_u_cnstr,1);
        u_S = zeros(num_u_cnstr,num_u);
        
        u_s_(1:num_u_cnstr/2) = ...
            [u_trans_cnstr*ones(3,1); u_rot_cnstr*ones(3,1)];
        u_s_(num_u_cnstr/2+1:end) = ...
            [u_trans_cnstr*ones(3,1); u_rot_cnstr*ones(3,1)];

        u_S(1:num_u_cnstr/2,:) = diag(ones(1,6));
        u_S(num_u_cnstr/2+1:end,:) = -diag(ones(1,6)); 
        U = Polyhedron('A',u_S,'B',u_s_);
        
        % Start algorithm to find Y_inf_beta approximation-----------------
        %------------------------------------------------------------------
                
        % Calculate Upsilon_k to get confidence ellipsoid
        Upsilon_ = zeros(k_prop,num_measure,num_measure);
        Xi_ = zeros(k_prop+1,12,12);
        Xi_(1,:,:) = Xi_0;
        Uncert_prop_ = zeros(k_prop,num_output,num_output);
        U_Uncert_prop_ = zeros(k_prop,num_u,num_u);
        
        for i = 1:k_prop
            Upsilon_(i,:,:) = H*squeeze(Xi_(i,:,:))*H'+V*P_NN*V';
            Xi_(i+1,:,:) = A*squeeze(Xi_(i,:,:))*A'+G*P_NN*G';
            Uncert_prop_(i,1:num_measure,1:num_measure) = Upsilon_(i,:,:);
            Uncert_prop_(i,num_measure+1:num_measure+12,num_measure+1:num_measure+12)...
                = Xi_(i,:,:);
            U_Uncert_prop_(i,:,:) = K*squeeze(Xi_(i,:,:))*K';
        end
        
        % Calculate Y_tilde_inf_beta 
        % Calculate U_tilde_inf_beta
        Y_tilde_inf_beta = Y;
        U_tilde_inf_beta = U;
        for i = 1:k_prop
            for ii = 1:num_cnstr
                ST = squeeze(S(ii,:));
                s_k(ii) = s_(ii)-sqrt(chi2inv(betaa,sum(ST~=0))).*...
                    sqrt(ST*squeeze(Uncert_prop_(i,:,:))*ST');
            end
            Y_k_beta = Polyhedron('A',S,'B',s_k);
            Y_tilde_inf_beta = intersect(Y_k_beta,Y_tilde_inf_beta);
            
            for ii = 1:num_u_cnstr
            u_ST = squeeze(u_S(ii,:));
            u_s_k(ii) = u_s_(ii)-sqrt(chi2inv(betaa,sum(u_ST~=0))).*...
                sqrt(u_ST*squeeze(U_Uncert_prop_(i,:,:))*u_ST');
            end
            U_k_beta = Polyhedron('A',u_S,'B',u_s_k);
            U_tilde_inf_beta = intersect(U_k_beta,U_tilde_inf_beta);
        end
        
        assignin('base','Y_tilde_inf_beta',Y_tilde_inf_beta);
        assignin('base','U',U_tilde_inf_beta);

        % Get rid of redundancies
        
        if k_prop ~= 0
            try
            Y_tilde_inf_beta = Y_tilde_inf_beta.minHRep();
            catch
                disp('    WARNING: Could not reduce redundancies. (first attempt)');
            end
            
            % Ensure not empty
            if isEmptySet(Y_tilde_inf_beta)
                disp('    WARNING: Y_tilde_inf_beta is empty.');
            end
            
            if isEmptySet(U_tilde_inf_beta)
                disp('    WARNING: U_tilde_inf_beta is empty.');
            end
            
        end
    
    end
        
end
    

    
    
    
    
    
    
    