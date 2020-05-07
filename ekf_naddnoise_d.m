function [xi_hat_t_t,Sig_t_t,y_guess,A_t,H_t,V_t] = ekf_naddnoise_d(t,y,u,xi,Sig,...
                                                        P_NN,Ffun,Hfun,HwrNfun,fp,f,tau)
    % Inputs and outputs
    % y             = Measured output (Expects new inputs)             15x1
    % xi            = [R V Th Omega_b]'  (Estimate) (Expects new inputs) 
    %                                                                  12x1
    % Sig           = Covariance of State Error (Expects new inputs)  
    %                                                                 12x12
    % R_t           = Covariance of Measurement Noise                 16x16
    % Ffun          = Function handle for Jacobian Evaluation of 
    %                 df/dxi. evaluates to A
    % Hfun          = Function handle for Jacobian Evaluation of dh/dxi
    %                 evaluates to H
    % HwrNfun       = Function handle for Jacobian Evaluation of dh/dN
    %                 evaluates to V
    % f             = focal length                                      1x1
    % tau           = sampling time                                     1x1
    % fp            = feature points [fp_1 ... fp_4]                    1x4 
    %                                          struct with fields x,y,z
    % xi_hat_t_t    = State estimate after this iteration              12x1
    % Sig_t_t       = Covariance of error estimate after this iteration
    %                                                                 12x12
    % A_t,H_t,V_t   = Linearizations. See control_alg for full description

    
    
    
    xi_hat_tm1_tm1 = xi;
    u_tm1 = u;
    Sig_tm1_tm1 = Sig;
    
    
    % W_t*Q_t_W_t' modeled as 
    necessary_process_noise = 0.000001.*eye(12);

    % Predict.
    xi_hat_t_tm1 = cntr_model_d(xi_hat_tm1_tm1,u_tm1,tau);
    
    A_t = Ffun(xi_hat_t_tm1(1),xi_hat_t_tm1(2),xi_hat_t_tm1(3),...
               xi_hat_t_tm1(8),xi_hat_t_tm1(9),xi_hat_t_tm1(10),...
               xi_hat_t_tm1(11),xi_hat_t_tm1(12));
           
    Sig_t_tm1 = A_t*Sig_tm1_tm1*A_t'+necessary_process_noise;
    
    H_t = Hfun(0,0,0,0,0,0,0,0,0,0,0,0,0,t,...
               xi_hat_t_tm1(1),xi_hat_t_tm1(2),xi_hat_t_tm1(3),...
               xi_hat_t_tm1(7),xi_hat_t_tm1(8),xi_hat_t_tm1(9));
    V_t = HwrNfun(0,0,0,0,0,0,0,0,0,0,0,0,0,t,...
               xi_hat_t_tm1(1),xi_hat_t_tm1(2),xi_hat_t_tm1(3),...
               xi_hat_t_tm1(7),xi_hat_t_tm1(8),xi_hat_t_tm1(9));
      
    % Update.
    L_t = Sig_t_tm1*H_t'/(H_t*Sig_t_tm1*H_t'+V_t*P_NN*V_t');
    y_guess = measure_model_d(t,xi_hat_t_tm1,fp,f,...
                               zeros(3,4),0,zeros(3,1));
    xi_hat_t_t = xi_hat_t_tm1+...
        L_t*(y-y_guess);
    Sig_t_t = (eye(12)-L_t*H_t)*Sig_t_tm1;
    
    
    % Check for errors
    if isnan(sum(A_t,'all'))
        disp('    ERROR: A_t has NaN');
    elseif isnan(sum(H_t,'all'))
        disp('    ERROR: H_t has NaN');
    elseif isnan(sum(V_t,'all'))
        disp('    ERROR: V_t has NaN');
    elseif isnan(sum(L_t,'all'))
        disp('    ERROR: L_t has NaN');
    end
   

end
