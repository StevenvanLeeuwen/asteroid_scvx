function [Y_tilde_inf_beta,U] = compute_constraints(Y_update,foview_subtract,...
                                A,G,H,V,...
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
    z_cnstr = 5.8;
    radius_cnstr = 20;
    u_trans_cnstr = 0.003;
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
        
        u_s_ = zeros(num_u_cnstr,1);
        u_S_T = zeros(num_u_cnstr,num_u);
        u_s_(1:num_u_cnstr/2) = ...
            [u_trans_cnstr*ones(3,1); u_rot_cnstr*ones(3,1)];
        u_s_(num_u_cnstr/2+1:end) = ...
            [u_trans_cnstr*ones(3,1); u_rot_cnstr*ones(3,1)];

        u_S_T(1:num_u_cnstr/2,:) = diag(ones(1,6));
        u_S_T(num_u_cnstr/2+1:end,:) = -diag(ones(1,6)); 
        U = Polyhedron('A',u_S_T,'B',u_s_);
        
        % Start algorithm to find Y_inf_beta approximation-----------------
        %------------------------------------------------------------------
                
        % Calculate Upsilon_k to get confidence ellipsoid
        Upsilon_ = zeros(k_prop,num_measure,num_measure);
        Xi_ = zeros(k_prop+1,12,12);
        Xi_(1,:,:) = Xi_0;
        Uncert_prop_ = zeros(k_prop,num_output,num_output);
        
        for i = 1:k_prop
            Upsilon_(i,:,:) = H*squeeze(Xi_(i,:,:))*H'+V*P_NN*V';
            Xi_(i+1,:,:) = A*squeeze(Xi_(i,:,:))*A'+G*P_NN*G';
            Uncert_prop_(i,1:num_measure,1:num_measure) = Upsilon_(i,:,:);
            Uncert_prop_(i,num_measure+1:num_measure+12,num_measure+1:num_measure+12)...
                = Xi_(i,:,:);
        end

        % Calculate Y_tilde_inf_beta 
        Y_tilde_inf_beta = Y;
        for i = 1:k_prop
            for ii = 1:num_cnstr
                ST = squeeze(S(ii,:));
                s_k(ii) = s_(ii)-sqrt(chi2inv(betaa,sum(ST~=0))).*...
                    sqrt(ST*squeeze(Uncert_prop_(i,:,:))*ST');
            end
            Y_k_beta = Polyhedron('A',S,'B',s_k);
            Y_tilde_inf_beta = intersect(Y_k_beta,Y_tilde_inf_beta);
        end
        
        assignin('base','Y_tilde_inf_beta',Y_tilde_inf_beta);

        % Get rid of redundancies
        A = 1;
        
        if k_prop ~= 0
            try
            [A,b] = noredund(Y_tilde_inf_beta.A(:,[5:12,num_measure+3]),...
                Y_tilde_inf_beta.b);
            catch
                disp('    WARNING: Could not reduce redundancies. (first attempt)');
            end
            
            A_noredund = zeros(size(A,1),num_output);
            A_noredund(:,[5:12,num_measure+3]) = A;
            b_noredund = b;
            Y_tilde_inf_beta = Polyhedron('A',A_noredund,'B',b_noredund);
            
            % Ensure not empty
            if isEmptySet(Y_tilde_inf_beta)
                disp('    WARNING: Y_tilde_inf_beta is empty.');
            end
            
        end
    
    end
        
end
    

    
    
    
    
    
    
    