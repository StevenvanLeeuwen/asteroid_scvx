

% THIS ONE ASSUMES A LINEAR SYSTEM ASTACK,CSTACK (not at 0), SO 
% NO REFERENCE TRACKING AND INIT COND IS NONZERO

function [ubar,xibar,Jval,feasible,output] = scvx_cnstr_slack_not_at_zero(...
    xistack,xirefstack,urefstack,Astack,B,Cstack,Y,U,N)

    %- - - - - - - - - - - - - - - - - - - - - - - - - -                        
    num_x = 12;
    num_u = 6;
    %- - - - - - - - - - - - - - - - - - - - - - - - - -  
    num_cnstr = size(Y.A,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constants
    R = 1*eye(num_u);                                    % control penalty
    R(1:3,1:3) = 5000*eye(3);
    R(3,3) = 5000;
    Q = 10*eye(num_x);                                 % state penalty
    Q(4:6,4:6) = 0.0001*eye(3);
    Q(10:12,10:12) = 0.0001*eye(3);
    mu = 10^9*eye(num_cnstr);                          % slack penalty
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Build the QP
    S = zeros(num_x*N,num_u*N);
    M = zeros(num_x*N,1);
    Cdiag = zeros(num_cnstr*N,num_x*N);
    W = zeros(num_cnstr*N,1);

    for cnt = 1:N 
        % Build the objective function components
        for k = 1:cnt
            if k == cnt
                S((k-1)*num_x+1:k*num_x,(k-1)*num_u+1:k*num_u) = B;  
            else
                S((cnt-1)*num_x+1:cnt*num_x,(k-1)*num_u+1:k*num_u) = ...
                    Astack((cnt-1)*num_x+1:cnt*num_x,:)*...
                    S((cnt-2)*num_x+1:(cnt-1)*num_x,(k-1)*num_u+1:k*num_u); 
            end
        end
        if cnt == 1
            M((cnt-1)*num_x+1:cnt*num_x,:) = Astack((cnt-1)*num_x+1:cnt*num_x,:)*xistack(1:num_x);
        else
            M((cnt-1)*num_x+1:cnt*num_x,:) = (Astack((cnt-1)*num_x+1:cnt*num_x,:)^cnt)*...
                xistack((cnt-1)*num_x+1:cnt*num_x);
        end

        Qbar((cnt-1)*num_x+1:cnt*num_x,(cnt-1)*num_x+1:cnt*num_x) = Q;
        Rbar((cnt-1)*num_u+1:cnt*num_u,(cnt-1)*num_u+1:cnt*num_u) = R;
        
        % Build the constraint components
        Cdiag((cnt-1)*num_cnstr+1:cnt*num_cnstr,(cnt-1)*num_x+1:cnt*num_x) =...
            Cstack((cnt-1)*num_cnstr+1:cnt*num_cnstr,:);
        W((cnt-1)*num_cnstr+1:cnt*num_cnstr) = Y.b;
    end
        G = Cdiag*S;
        W = W-(Cdiag*(M+xirefstack));

    % Add slack
    S = [S zeros(size(S,1),num_cnstr)];
    Rbar = [Rbar zeros(size(Rbar,1),num_cnstr);zeros(num_cnstr,size(Rbar,2)) mu];
    G_right = zeros(size(G,1),num_cnstr);
    for cnt = 1:N
        G_right((cnt-1)*num_cnstr+1:cnt*num_cnstr,:) = -eye(num_cnstr);
    end
    G = [G G_right];
    
    % Add control constraints
    G = [G;eye(N*num_u) zeros(N*num_u,num_cnstr);-eye(N*num_u) zeros(N*num_u,num_cnstr)];
    W_below1 = zeros(N*num_u,1);
    W_below2 = zeros(N*num_u,1);
    for cnt = 1:N
        W_below1((cnt-1)*num_u+1:cnt*num_u,:) = ...
            U.b(1:num_u)-urefstack((cnt-1)*num_u+1:cnt*num_u,:);
        W_below2((cnt-1)*num_u+1:cnt*num_u,:) = ...
            U.b(1:num_u)+urefstack((cnt-1)*num_u+1:cnt*num_u,:);
    end
    W = [W;W_below1;W_below2];
    
    % Solve
    H = 2*(S'*Qbar*S+Rbar);
    qT = 2*M'*Qbar*S;

    % quadprog(H,f,A,b,B,c,l,u,x0,options)
    options = optimset('Display','off');
    [ubar,Jval,feasible,output,~] = quadprog(H,qT',G,W,[],[],[],[],[],options);
    xibar = S*ubar+M+xirefstack;
    ubar = ubar(1:num_u*N)+urefstack;

end
    
        
    
      
