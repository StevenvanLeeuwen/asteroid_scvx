function [N_,N_c,N_b,P_NN] = measurement_noise_gen(Sigma,Cor)
    % Inputs and outputs
    % Sigma         = standard deviation of 
    %                  [N1...N4,N_c,N_b]'                              16x1
    % Cor           = Correlation coeff  matrix of 
    %                  [N1...N4,N_c,N_b]'                             16x16
    % N_,N_c,N_b    = noise samples
    % P_NN          = Covariance Matrix of [N1...N4,N_c,N_b]          16x16
    
    lenN = 16;
    s_mat = zeros(lenN,lenN);
    for k = 1:lenN
        for kk = 1:lenN
            s_mat(k,kk) = Sigma(k)*Sigma(kk);
        end
    end 
    P_NN = Cor.*s_mat;
    N = normrnd(0,1,[1,lenN])*chol(P_NN);
    
    N_ = [N(1:3)' N(4:6)' N(7:9)' N(10:12)'];
    N_c = N(13);
    N_b = N(14:16)';
    
end