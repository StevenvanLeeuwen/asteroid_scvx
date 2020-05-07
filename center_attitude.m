function ref_adj = center_attitude(ref)

    % Inputs and Outputs
    % ref                    = some reference state vector             12x1
    % ref_adj                = ref with euler angles recomputed        12x1
    
    % Center the ref euler angles to the average fp location
    
    % Constants
    omega=2*pi/(5.27*3600);        %rotation rate (rads/sec)

    t = evalin('base','t');
    ref_adj_prev = evalin('base','ref_adj');
    fp = evalin('base','fp');
    rnd1 = evalin('base','rnd1');
    rnd2 = evalin('base','rnd2');
    
    pi_tol = 0.1;
    
    avg_fp = mean([[fp(1).x;fp(2).x;fp(3).x;fp(4).x],...
         [fp(1).y;fp(2).y;fp(3).y;fp(4).y],...
         [fp(1).z;fp(2).z;fp(3).z;fp(4).z]])';
    
    Ria = [cos(omega*t) -sin(omega*t) 0;...
           sin(omega*t) cos(omega*t) 0;...
           0 0 1];
    
    % Compute adjusted reference
    ref(7:9) = euler_angle_calc(Ria*(ref(1:3)-avg_fp),rnd1,rnd2);
      
    for k = 7:9
        if abs(ref(k)-ref_adj_prev(k)) > pi_tol
            if abs(ref(k)+pi-ref_adj_prev(k)) <= pi_tol
                ref(k) = ref(k)+pi;
            elseif abs(ref(k)-pi-ref_adj_prev(k)) <= pi_tol
                ref(k) = ref(k)-pi;
            elseif abs(ref(k)+2*pi-ref_adj_prev(k)) <= pi_tol
                ref(k) = ref(k)+2*pi;
            elseif abs(ref(k)-2*pi-ref_adj_prev(k)) <= pi_tol
                ref(k) = ref(k)-2*pi;
            end
        end
    end
    
    ref_adj = ref;
    
end 

%     omega = evalin('base','omega');
%     t = evalin('base','t');
%     % This is the convention of how our Euler Angles are defined.
%         % Below must match that in init_focal_length.m
%         Ria = [cos(omega*t) -sin(omega*t) 0;...
%                sin(omega*t) cos(omega*t) 0;...
%                0 0 1];
%         R1i = [cos(ref_adj(9)) -sin(ref_adj(9)) 0;...
%                sin(ref_adj(9)) cos(ref_adj(9)) 0;...
%                0 0 1];
%         R21 = [cos(ref_adj(8)) 0 sin(ref_adj(8));...
%                0 1 0;...
%                -sin(ref_adj(8)) 0 cos(ref_adj(8))];
%         Rb2 = [1 0 0;...
%                0 cos(ref_adj(7)) -sin(ref_adj(7));...
%                0 sin(ref_adj(7)) cos(ref_adj(7))];
% 
%         Rba = Rb2*R21*R1i*Ria;
% 
%         % Want the transformation from eul to rotm to be accurate according to
%         % tol
%         check = Rba*(xi_0(1:3)-avg_fp)





