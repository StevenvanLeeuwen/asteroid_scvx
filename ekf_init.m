function Sig = ekf_init()
     Sig = 0.00000002*eye(12);
%      Sig(1,1) = 0.7;
%      Sig(2,2) = 0.6;
%      Sig(3,3) = 0.02;
%      Sig(4,4) = 0.4;
%      Sig(5,5) = 0.002;
%      Sig(6,6) = 0.002;
%      Sig(7,7) = 0.3;
%      Sig(8,8) = 0.4;
%      Sig(9,9) = 0.3;
%      Sig(10,10) = 0.3;
%      Sig(11,11) = 0.4;
%      Sig(12,12) = 0.3;
end
