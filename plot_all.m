%% Plot
disp('Plotting...');

disp('Making 3D asteroid...');
tdfread('eros001708_WD.tab');       %Reads the shape data
[vert,faces]=Eros_setup(vf,x,y,z);  %Sorts the data into vertices and faces

% Plot 3D track, as well as feature points and starting point.
figure
plot3(xi_record(1,1),xi_record(1,2),xi_record(1,3),'b*','MarkerSize',20,'LineWidth',2)
hold on
plot3(xi_record(end,1),xi_record(end,2),xi_record(end,3),'bo','MarkerSize',20,'LineWidth',2)
plot3(fp(1).x,fp(1).y,fp(1).z,'.r','MarkerSize',50,'LineWidth',2)
plot3(fp(2).x,fp(2).y,fp(2).z,'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot3(fp(3).x,fp(3).y,fp(3).z,'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot3(fp(4).x,fp(4).y,fp(4).z,'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot3(xi_record(:,1),xi_record(:,2),xi_record(:,3),'b','LineWidth',2)
plot3(xi_nonvis_record(:,1),xi_nonvis_record(:,2),xi_nonvis_record(:,3),...
    'k--','LineWidth',2)
plot3(xi_nonvis_record(end,1),xi_nonvis_record(end,2),xi_nonvis_record(end,3),...
    'ko','MarkerSize',20,'LineWidth',2)
h=patch('vertices',vert,'faces',faces,'facecolor',[0.95,0.95,0.95]);
xlabel('$x$ (km)','Interpreter','latex')
ylabel('$y$ (km)','Interpreter','latex')
zlabel('$z$ (km)','Interpreter','latex')
legend('Starting Point','Goal','Feature Point','Estimated Position',...
    'Actual Position','Actual End Position');
title('Track of Position','Interpreter','latex');
axis equal
grid on
%

% Plot histories r
figure
subplot(3,1,1)
plot(tspan,xi_record(:,1),'b','LineWidth',2)
hold on
plot(tspan,xi_nonvis_record(:,1),'k','LineWidth',2)
plot(tspan,ref(1,:),'k--','LineWidth',2)
hold off
legend('estimate','actual','ref','location','northeast')
ylabel('$x$ (km)','Interpreter','latex')
title('History of $\mathbf{r}$','Interpreter','latex')
grid on

subplot(3,1,2)
plot(tspan,xi_record(:,2),'b','LineWidth',2)
hold on
plot(tspan,xi_nonvis_record(:,2),'k','LineWidth',2)
plot(tspan,ref(2,:),'k--','LineWidth',2)
hold off
legend('estimate','actual','ref','location','northeast')
ylabel('$y$ (km)','Interpreter','latex')
grid on

subplot(3,1,3)
plot(tspan,xi_record(:,3),'b','LineWidth',2)
hold on
plot(tspan,xi_nonvis_record(:,3),'k','LineWidth',2)
plot(tspan,ref(3,:),'k--','LineWidth',2)
plot(tspan,z_cnstr.*ones(1,length(tspan)),'r--','LineWidth',2)
hold off
legend('estimate','actual','ref','constraint')
xlabel('time (s)','Interpreter','latex');
ylabel('$z$ (km)','Interpreter','latex');
grid on
%
 
% Plot histories Th
figure
subplot(3,1,1)
plot(tspan,xi_record(:,7),'b','LineWidth',2)
hold on
plot(tspan,xi_nonvis_record(:,7),'k','LineWidth',2)
plot(tspan,ref_adj_record(7,:),'k--','LineWidth',2)
%plot(tspan,ref_adj_record(7,:),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
hold off
legend('estimate','actual','ref')
ylabel('$\phi$ (rad)','Interpreter','latex')
title('History of $\mathbf{\Theta}$','Interpreter','latex')
grid on

subplot(3,1,2)
plot(tspan,xi_record(:,8),'b','MarkerSize',20,'LineWidth',2)
hold on
plot(tspan,xi_nonvis_record(:,8),'k','LineWidth',2)
plot(tspan,ref_adj_record(8,:),'k--','LineWidth',2)
%plot(tspan,ref_adj_record(8,:),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
hold off
legend('estimate','actual','ref')
ylabel('$\theta$ (rad)','Interpreter','latex')
grid on

subplot(3,1,3)
plot(tspan,xi_record(:,9),'b','MarkerSize',20,'LineWidth',2)
hold on
plot(tspan,xi_nonvis_record(:,9),'k','LineWidth',2)
plot(tspan,ref_adj_record(9,:),'k--','LineWidth',2)
%plot(tspan,ref_adj_record(9,:),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
hold off
legend('estimate','actual','ref')
xlabel('time (s)');
ylabel('$\psi$ (rad)','Interpreter','latex');
grid on
%

% Plot histories Angular Velocity
figure
subplot(3,1,1)
plot(tspan,xi_record(:,10),'b','LineWidth',2)
hold on
plot(tspan,xi_nonvis_record(:,10),'k','LineWidth',2)
plot(tspan,ref(10,:),'k--','LineWidth',2)
plot(tspan,ref_adj_record(10,:),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
hold off
legend('estimate','actual','ref','ref adj')
ylabel('$\omega_{b,1}$ (rad/s)','Interpreter','latex')
title('History of $\mathbf{\omega_b}$','Interpreter','latex')
grid on

subplot(3,1,2)
plot(tspan,xi_record(:,11),'b','MarkerSize',20,'LineWidth',2)
hold on
plot(tspan,xi_nonvis_record(:,11),'k','LineWidth',2)
plot(tspan,ref(11,:),'k--','LineWidth',2)
plot(tspan,ref_adj_record(11,:),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
hold off
legend('estimate','actual','ref','ref adj')
ylabel('$\omega_{b,2}$ (rad/s)','Interpreter','latex')
grid on

subplot(3,1,3)
plot(tspan,xi_record(:,12),'b','MarkerSize',20,'LineWidth',2)
hold on
plot(tspan,xi_nonvis_record(:,12),'k','LineWidth',2)
plot(tspan,ref(12,:),'k--','LineWidth',2)
plot(tspan,ref_adj_record(12,:),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
hold off
legend('estimate','actual','ref','ref adj')
xlabel('time (s)');
ylabel('$\omega_{b,3}$ (rad/s)','Interpreter','latex');
grid on
%

% Plot histories u
utrans_min = -u_trans_cnstr;
utrans_max = u_trans_cnstr;

urot_min = -u_rot_cnstr;
urot_max = u_rot_cnstr;

figure
subplot(2,1,1)
plot(tspan,u_record(:,1),'b','LineWidth',2);
hold on
plot(tspan,u_record(:,2),'c','LineWidth',2);
plot(tspan,u_record(:,3),'m','LineWidth',2);
plot(tspan,utrans_min*ones(1,length(tspan)),'r--','LineWidth',2)
plot(tspan,utrans_max*ones(1,length(tspan)),'r--','LineWidth',2)
hold off
legend('$u_x$','$u_y$','$u_z$','constraint','Interpreter','latex')
ylabel('(km/s)','Interpreter','latex')
title('History of $\mathbf{u}$','Interpreter','latex')
grid on

subplot(2,1,2)
plot(tspan,u_record(:,4),'b','LineWidth',2);
hold on
plot(tspan,u_record(:,5),'c','LineWidth',2);
plot(tspan,u_record(:,6),'m','LineWidth',2);
plot(tspan,urot_min*ones(1,length(tspan)),'r--','LineWidth',2)
plot(tspan,urot_max*ones(1,length(tspan)),'r--','LineWidth',2)
hold off
legend('$M_1$','$M_2$','$M_3$','constraint','Interpreter','latex')
xlabel('time (s)')
ylabel('$(\frac{kgm^2}{s^2})$','Interpreter','latex')
title('History of $\mathbf{M}$','Interpreter','latex')
grid on
%

% Plot histories y that are camera
figure
plot(y_record(:,5),y_record(:,6),'r','LineWidth',2)
hold on
plot(y_record(:,7),y_record(:,8),'r','LineWidth',2,'HandleVisibility','off')
plot(y_record(:,9),y_record(:,10),'r','LineWidth',2,'HandleVisibility','off')
plot(y_record(:,11),y_record(:,12),'r','LineWidth',2,'HandleVisibility','off')
plot(y_record(1,5),y_record(1,6),'.r','MarkerSize',50,'LineWidth',2)
plot(y_record(1,7),y_record(1,8),'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot(y_record(1,9),y_record(1,10),'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot(y_record(1,11),y_record(1,12),'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot([radius_cnstr -radius_cnstr -radius_cnstr radius_cnstr radius_cnstr],...
    [radius_cnstr radius_cnstr -radius_cnstr -radius_cnstr radius_cnstr],...
    'r--','LineWidth',2);
hold off
legend('Trace of feature points','Starting Point',...
    'Constraints','location','Southeast')
xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')
axis equal
title('History of $\mathbf{c_1}...\mathbf{c_4}$','Interpreter','latex')
grid on

% Plot histories y that are camera
figure
plot(y_fromref_record(:,5),y_fromref_record(:,6),'r','LineWidth',2)
hold on
plot(y_fromref_record(:,7),y_fromref_record(:,8),'r','LineWidth',2,'HandleVisibility','off')
plot(y_fromref_record(:,9),y_fromref_record(:,10),'r','LineWidth',2,'HandleVisibility','off')
plot(y_fromref_record(:,11),y_fromref_record(:,12),'r','LineWidth',2,'HandleVisibility','off')
plot(y_fromref_record(1,5),y_fromref_record(1,6),'.r','MarkerSize',50,'LineWidth',2)
plot(y_fromref_record(1,7),y_fromref_record(1,8),'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot(y_fromref_record(1,9),y_fromref_record(1,10),'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot(y_fromref_record(1,11),y_fromref_record(1,12),'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot([radius_cnstr -radius_cnstr -radius_cnstr radius_cnstr radius_cnstr],...
    [radius_cnstr radius_cnstr -radius_cnstr -radius_cnstr radius_cnstr],...
    'r--','LineWidth',2);
hold off
legend('Trace of feature points','Starting Point',...
    'Constraints','location','Southeast')
xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')
axis equal
title('History of $\mathbf{c_1}...\mathbf{c_4} \ OF \ REF \ Note \ depends \ on \ rnd1,rnd2 $','Interpreter','latex')
grid on

% Plot state errors

figure
subplot(4,3,1);
plot(tspan,abs(xi_record(:,1)-xi_nonvis_record(:,1)),'b','LineWidth',2);
ylabel('$x$ (km)','Interpreter','latex')
sgtitle('Estimate Error')
grid on

subplot(4,3,2);
plot(tspan,abs(xi_record(:,2)-xi_nonvis_record(:,2)),'b','LineWidth',2);
ylabel('$y$ (km)','Interpreter','latex')
grid on

subplot(4,3,3);
plot(tspan,abs(xi_record(:,3)-xi_nonvis_record(:,3)),'b','LineWidth',2);
ylabel('$z$ (km)','Interpreter','latex')
grid on

subplot(4,3,4);
plot(tspan,abs(xi_record(:,4)-xi_nonvis_record(:,4)),'b','LineWidth',2);
ylabel('$\dot{x}$ (km/s)','Interpreter','latex')
grid on

subplot(4,3,5);
plot(tspan,abs(xi_record(:,5)-xi_nonvis_record(:,5)),'b','LineWidth',2);
ylabel('$\dot{y}$ (km/s)','Interpreter','latex')
grid on

subplot(4,3,6);
plot(tspan,abs(xi_record(:,6)-xi_nonvis_record(:,6)),'b','LineWidth',2);
ylabel('$\dot{z}$ (km/s)','Interpreter','latex')
grid on

subplot(4,3,7);
plot(tspan,abs(xi_record(:,7)-xi_nonvis_record(:,7)),'b','LineWidth',2);
ylabel('$\phi$ (rad)','Interpreter','latex')
grid on

subplot(4,3,8);
plot(tspan,abs(xi_record(:,8)-xi_nonvis_record(:,8)),'b','LineWidth',2);
ylabel('$\theta$ (rad)','Interpreter','latex')
grid on

subplot(4,3,9);
plot(tspan,abs(xi_record(:,9)-xi_nonvis_record(:,9)),'b','LineWidth',2);
ylabel('$\psi$ (rad)','Interpreter','latex')
grid on

subplot(4,3,10);
plot(tspan,abs(xi_record(:,10)-xi_nonvis_record(:,10)),'b','LineWidth',2);
xlabel('time (s)')
ylabel('$\omega_x$ (rad/s)','Interpreter','latex')
grid on

subplot(4,3,11);
plot(tspan,abs(xi_record(:,11)-xi_nonvis_record(:,11)),'b','LineWidth',2);
xlabel('time (s)')
ylabel('$\omega_y$ (rad/s)','Interpreter','latex')
grid on

subplot(4,3,12);
plot(tspan,abs(xi_record(:,12)-xi_nonvis_record(:,12)),'b','LineWidth',2);
xlabel('time (s)')
ylabel('$\omega_z$ (rad/s)','Interpreter','latex')
grid on

% Plot covariance errors between state measurements and prediction

figure
subplot(4,3,1);
plot(tspan,Sig_record(:,1,1),'b','LineWidth',2);
ylabel('$x$ (km)','Interpreter','latex')
sgtitle('Covariance Error EKF')
grid on

subplot(4,3,2);
plot(tspan,Sig_record(:,2,2),'b','LineWidth',2);
ylabel('$y$ (km)','Interpreter','latex')
grid on

subplot(4,3,3);
plot(tspan,Sig_record(:,3,3),'b','LineWidth',2);
ylabel('$z$ (km)','Interpreter','latex')
grid on

subplot(4,3,4);
plot(tspan,Sig_record(:,4,4),'b','LineWidth',2);
ylabel('$\dot{x}$ (km/s)','Interpreter','latex')
grid on

subplot(4,3,5);
plot(tspan,Sig_record(:,5,5),'b','LineWidth',2);
ylabel('$\dot{y}$ (km/s)','Interpreter','latex')
grid on

subplot(4,3,6);
plot(tspan,Sig_record(:,6,6),'b','LineWidth',2);
ylabel('$\dot{z}$ (km/s)','Interpreter','latex')
grid on

subplot(4,3,7);
plot(tspan,Sig_record(:,7,7),'b','LineWidth',2);
ylabel('$\phi$ (rad)','Interpreter','latex')
grid on

subplot(4,3,8);
plot(tspan,Sig_record(:,8,8),'b','LineWidth',2);
ylabel('$\theta$ (rad)','Interpreter','latex')
grid on

subplot(4,3,9);
plot(tspan,Sig_record(:,9,9),'b','LineWidth',2);
ylabel('$\psi$ (rad)','Interpreter','latex')
grid on

subplot(4,3,10);
plot(tspan,Sig_record(:,10,10),'b','LineWidth',2);
xlabel('time (s)')
ylabel('$\omega_x$ (rad/s)','Interpreter','latex')
grid on

subplot(4,3,11);
plot(tspan,Sig_record(:,11,11),'b','LineWidth',2);
xlabel('time (s)')
ylabel('$\omega_y$ (rad/s)','Interpreter','latex')
grid on

subplot(4,3,12);
plot(tspan,Sig_record(:,12,12),'b','LineWidth',2);
xlabel('time (s)')
ylabel('$\omega_z$ (rad/s)','Interpreter','latex')
grid on

% Plot measurement discrepency

figure
subplot(4,4,1);
plot(tspan,abs(y_record(:,1)-y_guess_record(:,1)),'b','LineWidth',2);
sgtitle('Measurement Discrepency y and prediction')
grid on

subplot(4,4,2);
plot(tspan,abs(y_record(:,2)-y_guess_record(:,2)),'b','LineWidth',2);
grid on

subplot(4,4,3);
plot(tspan,abs(y_record(:,3)-y_guess_record(:,3)),'b','LineWidth',2);
grid on

subplot(4,4,4);
plot(tspan,abs(y_record(:,4)-y_guess_record(:,4)),'b','LineWidth',2);
grid on

subplot(4,4,5);
plot(tspan,abs(y_record(:,5)-y_guess_record(:,5)),'b','LineWidth',2);
grid on

subplot(4,4,6);
plot(tspan,abs(y_record(:,6)-y_guess_record(:,6)),'b','LineWidth',2);
grid on

subplot(4,4,7);
plot(tspan,abs(y_record(:,7)-y_guess_record(:,7)),'b','LineWidth',2);
grid on

subplot(4,4,8);
plot(tspan,abs(y_record(:,8)-y_guess_record(:,8)),'b','LineWidth',2);
grid on

subplot(4,4,9);
plot(tspan,abs(y_record(:,9)-y_guess_record(:,9)),'b','LineWidth',2);
grid on

subplot(4,4,10);
plot(tspan,abs(y_record(:,10)-y_guess_record(:,10)),'b','LineWidth',2);
grid on

subplot(4,4,11);
plot(tspan,abs(y_record(:,11)-y_guess_record(:,11)),'b','LineWidth',2);
grid on

subplot(4,4,12);
plot(tspan,abs(y_record(:,12)-y_guess_record(:,12)),'b','LineWidth',2);
grid on

subplot(4,4,13);
plot(tspan,abs(y_record(:,13)-y_guess_record(:,13)),'b','LineWidth',2);
grid on

subplot(4,4,14);
plot(tspan,abs(y_record(:,14)-y_guess_record(:,14)),'b','LineWidth',2);
grid on

subplot(4,4,15);
plot(tspan,abs(y_record(:,15)-y_guess_record(:,15)),'b','LineWidth',2);
grid on

% % Plot measurement discrepency
% 
% figure
% subplot(4,4,1);
% plot(tspan,abs(y_record(:,1)-y_fromref_record(:,1)),'b','LineWidth',2);
% sgtitle('Measurement Discrepency y prediction and non noise case')
% grid on
% 
% subplot(4,4,2);
% plot(tspan,abs(y_guess_record(:,2)-y_fromref_record(:,2)),'b','LineWidth',2);
% grid on
% 
% subplot(4,4,3);
% plot(tspan,abs(y_guess_record(:,3)-y_fromref_record(:,3)),'b','LineWidth',2);
% grid on
% 
% subplot(4,4,4);
% plot(tspan,abs(y_guess_record(:,4)-y_fromref_record(:,4)),'b','LineWidth',2);
% grid on
% 
% subplot(4,4,5);
% plot(tspan,abs(y_guess_record(:,5)-y_fromref_record(:,5)),'b','LineWidth',2);
% grid on
% 
% subplot(4,4,6);
% plot(tspan,abs(y_guess_record(:,6)-y_fromref_record(:,6)),'b','LineWidth',2);
% grid on
% 
% subplot(4,4,7);
% plot(tspan,abs(y_guess_record(:,7)-y_fromref_record(:,7)),'b','LineWidth',2);
% grid on
% 
% subplot(4,4,8);
% plot(tspan,abs(y_guess_record(:,8)-y_fromref_record(:,8)),'b','LineWidth',2);
% grid on
% 
% subplot(4,4,9);
% plot(tspan,abs(y_guess_record(:,9)-y_fromref_record(:,9)),'b','LineWidth',2);
% grid on
% 
% subplot(4,4,10);
% plot(tspan,abs(y_guess_record(:,10)-y_fromref_record(:,10)),'b','LineWidth',2);
% grid on
% 
% subplot(4,4,11);
% plot(tspan,abs(y_guess_record(:,11)-y_fromref_record(:,11)),'b','LineWidth',2);
% grid on
% 
% subplot(4,4,12);
% plot(tspan,abs(y_guess_record(:,12)-y_fromref_record(:,12)),'b','LineWidth',2);
% grid on
% 
% subplot(4,4,13);
% plot(tspan,abs(y_guess_record(:,13)-y_fromref_record(:,13)),'b','LineWidth',2);
% grid on
% 
% subplot(4,4,14);
% plot(tspan,abs(y_guess_record(:,14)-y_fromref_record(:,14)),'b','LineWidth',2);
% grid on
% 
% subplot(4,4,15);
% plot(tspan,abs(y_guess_record(:,15)-y_fromref_record(:,15)),'b','LineWidth',2);
% grid on

% Err dist vs camera 
% err_dist = sqrt((xi_record(:,1)-xi_nonvis_record(:,1)).^2+...
% (xi_record(:,2)-xi_nonvis_record(:,2)).^2+(xi_record(:,3)-xi_nonvis_record(:,3)).^2);
% 
% cam_dist1 = sqrt(y_record(:,5).^2+y_record(:,6).^2)-...
%     (sqrt(y_record(1,5).^2+y_record(1,6).^2));
% 
% cam_dist2 = sqrt(y_record(:,7).^2+y_record(:,8).^2)-...
%     (sqrt(y_record(1,7).^2+y_record(1,8).^2));
% 
% cam_dist3 = sqrt(y_record(:,9).^2+y_record(:,10).^2)-...
%     (sqrt(y_record(1,9).^2+y_record(1,10).^2));
% 
% cam_dist4 = sqrt(y_record(:,11).^2+y_record(:,12).^2)-...
%      (sqrt(y_record(1,11).^2+y_record(1,12).^2));
% 
% cam_dist = max([cam_dist1 cam_dist2 cam_dist3 cam_dist4],[],2);

% figure
% plot(cam_dist,err_dist,'bo-','LineWidth',2);
% hold on
% plot(cam_dist(1),err_dist(1),'b*','LineWidth',2,'MarkerSize',20);
