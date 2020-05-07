%% Plot 2 segments

% In format
% ref1
% ref2
% tspan
% u_record1
% u_record2
% xi_nonvis_record1
% xi_nonvis_record2
% xi_record1
% xi_record2
% y_record1
% y_record2


disp('Plotting...');

disp('Making 3D asteroid...');
tdfread('eros001708_WD.tab');       %Reads the shape data
[vert,faces]=Eros_setup(vf,x,y,z);  %Sorts the data into vertices and faces
% Plot 3D track, as well as feature points and starting point.
%%
figure
subplot(2,1,1)
plot3(xi_record1(1,1),xi_record1(1,2),xi_record1(1,3),'b*','MarkerSize',20,'LineWidth',2)
hold on
plot3(xi_record1(end,1),xi_record1(end,2),xi_record1(end,3),'bo','MarkerSize',20,'LineWidth',2)
plot3(xi_record1(:,1),xi_record1(:,2),xi_record1(:,3),'b','LineWidth',2)
plot3(xi_nonvis_record1(:,1),xi_nonvis_record1(:,2),xi_nonvis_record1(:,3),...
    'k--','LineWidth',2)
plot3(xi_nonvis_record1(end,1),xi_nonvis_record1(end,2),xi_nonvis_record1(end,3),...
    'ko','MarkerSize',20,'LineWidth',2)
plot3(xi_record1(1,1),xi_record1(1,2),xi_record1(1,3),'b*','MarkerSize',20,'LineWidth',2)
load fp_close1.mat
plot3(fp(1).x,fp(1).y,fp(1).z,'.r','MarkerSize',50,'LineWidth',2)
plot3(fp(2).x,fp(2).y,fp(2).z,'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot3(fp(3).x,fp(3).y,fp(3).z,'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot3(fp(4).x,fp(4).y,fp(4).z,'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
h=patch('vertices',vert,'faces',faces,'facecolor',[0.95,0.95,0.95]);
hold off
xlabel('$x$ (km)','Interpreter','latex')
ylabel('$y$ (km)','Interpreter','latex')
zlabel('$z$ (km)','Interpreter','latex')
legend('Starting Point','Goal','Estimated Position',...
    'Actual Position','Actual End Position');
sgtitle('Track of Position','Interpreter','latex','FontSize',20);
axis equal
grid on
subplot(2,1,2)
plot3(xi_record2(1,1),xi_record2(1,2),xi_record2(1,3),'b*','MarkerSize',20,'LineWidth',2)
hold on
plot3(xi_record2(end,1),xi_record2(end,2),xi_record2(end,3),'bo','MarkerSize',20,'LineWidth',2)
plot3(xi_record2(:,1),xi_record2(:,2),xi_record2(:,3),'b','LineWidth',2)
plot3(xi_nonvis_record2(:,1),xi_nonvis_record2(:,2),xi_nonvis_record2(:,3),...
    'k--','LineWidth',2)
plot3(xi_nonvis_record2(end,1),xi_nonvis_record2(end,2),xi_nonvis_record2(end,3),...
    'ko','MarkerSize',20,'LineWidth',2)
plot3(xi_record2(1,1),xi_record2(1,2),xi_record2(1,3),'b*','MarkerSize',20,'LineWidth',2)
load fp_close.mat
plot3(fp(1).x,fp(1).y,fp(1).z,'.r','MarkerSize',50,'LineWidth',2)
plot3(fp(2).x,fp(2).y,fp(2).z,'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot3(fp(3).x,fp(3).y,fp(3).z,'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot3(fp(4).x,fp(4).y,fp(4).z,'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
h=patch('vertices',vert,'faces',faces,'facecolor',[0.95,0.95,0.95]);
hold off
xlabel('$x$ (km)','Interpreter','latex')
ylabel('$y$ (km)','Interpreter','latex')
zlabel('$z$ (km)','Interpreter','latex')
% legend('Starting Point','Goal','Estimated Position',...
%     'Actual Position','Actual End Position');
axis equal
grid on
%%
% Plot histories r
% figure
% subplot(3,2,1)
% plot(tspan1,xi_record1(:,1),'b','LineWidth',2)
% hold on
% plot(tspan1,xi_nonvis_record1(:,1),'k','LineWidth',2)
% plot(tspan1,ref1(1,:),'k--','LineWidth',2)
% hold off
% legend('estimate','actual','ref','location','northwest')
% ylabel('$x$ (km)','Interpreter','latex')
% sgtitle('History of $\mathbf{r}$','Interpreter','latex')
% title('Leg 1');
% grid on
% subplot(3,2,3)
% plot(tspan1,xi_record1(:,2),'b','LineWidth',2)
% hold on
% plot(tspan1,xi_nonvis_record1(:,2),'k','LineWidth',2)
% plot(tspan1,ref1(2,:),'k--','LineWidth',2)
% hold off
% legend('estimate','actual','ref','location','northeast')
% ylabel('$y$ (km)','Interpreter','latex')
% grid on
subplot(2,1,1)
plot(tspan1,xi_record1(:,3),'b','LineWidth',2)
hold on
plot(tspan1,xi_nonvis_record1(:,3),'k','LineWidth',2)
plot(tspan1,ref1(3,:),'k--','LineWidth',2)
plot(tspan1,z_cnstr1.*ones(1,length(tspan1)),'r--','LineWidth',2)
hold off
legend('estimate','actual','ref','constraint','location','northeast')
xlabel('time (s)','Interpreter','latex')
ylabel('$z$ (km)','Interpreter','latex')
grid on

% subplot(3,2,2)
% plot(tspan2,xi_record2(:,1),'b','LineWidth',2)
% hold on
% plot(tspan2,xi_nonvis_record2(:,1),'k','LineWidth',2)
% plot(tspan2,ref2(1,:),'k--','LineWidth',2)
% hold off
% legend('estimate','actual','ref','location','southeast')
% ylabel('$x$ (km)','Interpreter','latex')
% title('Leg 2');
% grid on
% subplot(3,2,4)
% plot(tspan2,xi_record2(:,2),'b','LineWidth',2)
% hold on
% plot(tspan2,xi_nonvis_record2(:,2),'k','LineWidth',2)
% plot(tspan2,ref2(2,:),'k--','LineWidth',2)
% hold off
% legend('estimate','actual','ref','location','northeast')
% ylabel('$y$ (km)','Interpreter','latex')
% grid on
subplot(2,1,2)
plot(tspan2,xi_record2(:,3),'b','LineWidth',2)
hold on
plot(tspan2,xi_nonvis_record2(:,3),'k','LineWidth',2)
plot(tspan2,ref2(3,:),'k--','LineWidth',2)
plot(tspan2,z_cnstr2.*ones(1,length(tspan2)),'r--','LineWidth',2)
hold off
%legend('estimate','actual','ref','constraint','location','northeast')
xlabel('time (s)','interpreter','latex')
ylabel('$z$ (km)','Interpreter','latex')
sgtitle('History of $z$','Interpreter','latex','FontSize',20)
grid on
%
%%
% Plot histories Th
figure
subplot(3,2,1)
plot(tspan1,xi_record1(:,7),'b','LineWidth',2)
hold on
plot(tspan1,xi_nonvis_record1(:,7),'k','LineWidth',2)
plot(tspan1,ref_adj_record1(7,:),'k--','LineWidth',2)
hold off
legend('estimate','actual','ref','location','northeast')
ylabel('$\phi$ (rad)','Interpreter','latex')
sgtitle('History of $\mathbf{\Theta}$','Interpreter','latex')
title('Leg 1')
grid on
subplot(3,2,3)
plot(tspan1,xi_record1(:,8),'b','LineWidth',2)
hold on
plot(tspan1,xi_nonvis_record1(:,8),'k','LineWidth',2)
plot(tspan1,ref_adj_record1(8,:),'k--','LineWidth',2)
hold off
legend('estimate','actual','ref','location','northwest')
ylabel('$\theta$ (rad)','Interpreter','latex')
grid on
subplot(3,2,5)
plot(tspan1,xi_record1(:,9),'b','LineWidth',2)
hold on
plot(tspan1,xi_nonvis_record1(:,9),'k','LineWidth',2)
plot(tspan1,ref_adj_record1(9,:),'k--','LineWidth',2)
hold off
legend('estimate','actual','ref','location','northwest')
xlabel('time (s)')
ylabel('$\psi$ (rad)','Interpreter','latex');
grid on

subplot(3,2,2)
plot(tspan2,xi_record2(:,7),'b','LineWidth',2)
hold on
plot(tspan2,xi_nonvis_record2(:,7),'k','LineWidth',2)
plot(tspan2,ref_adj_record2(7,:),'k--','LineWidth',2)
hold off
legend('estimate','actual','ref','location','southeast')
ylabel('$\phi$ (rad)','Interpreter','latex')
title('Leg 2')
grid on
subplot(3,2,4)
plot(tspan2,xi_record2(:,8),'b','LineWidth',2)
hold on
plot(tspan2,xi_nonvis_record2(:,8),'k','LineWidth',2)
plot(tspan2,ref_adj_record2(8,:),'k--','LineWidth',2)
hold off
legend('estimate','actual','ref','location','northeast')
ylabel('$\theta$ (rad)','Interpreter','latex')
grid on
subplot(3,2,6)
plot(tspan2,xi_record2(:,9),'b','LineWidth',2)
hold on
plot(tspan2,xi_nonvis_record2(:,9),'k','LineWidth',2)
plot(tspan2,ref_adj_record2(9,:),'k--','LineWidth',2)
hold off
legend('estimate','actual','ref','location','southeast')
xlabel('time (s)')
ylabel('$\psi$ (rad)','Interpreter','latex');
grid on
%
%%

% Plot histories u
utrans_min1 = -utrans_max1;
urot_min1 = -urot_max1;
utrans_min2 = -utrans_max2;
urot_min2 = -urot_max2;

figure
subplot(2,2,1)
plot(tspan1,u_record1(:,1),'b','LineWidth',2);
hold on
plot(tspan1,u_record1(:,2),'c','LineWidth',2);
plot(tspan1,u_record1(:,3),'m','LineWidth',2);
plot(tspan1,utrans_min1*ones(1,length(tspan1)),'r--','LineWidth',2)
plot(tspan1,utrans_max1*ones(1,length(tspan1)),'r--','LineWidth',2)
hold off
ylim([utrans_min1 -utrans_min1])
legend('$u_x$','$u_y$','$u_z$','constraint','Interpreter','latex')
ylabel('(km/s)','Interpreter','latex')
sgtitle('History of $\mathbf{u}$','Interpreter','latex')
title('Leg 1')
grid on

subplot(2,2,3)
plot(tspan1,u_record1(:,4),'b','LineWidth',2);
hold on
plot(tspan1,u_record1(:,5),'c','LineWidth',2);
plot(tspan1,u_record1(:,6),'m','LineWidth',2);
plot(tspan1,urot_min1*ones(1,length(tspan1)),'r--','LineWidth',2)
plot(tspan1,urot_max1*ones(1,length(tspan1)),'r--','LineWidth',2)
hold off
ylim([urot_min1 -urot_min1])
legend('$M_1$','$M_2$','$M_3$','constraint','Interpreter','latex')
xlabel('time (s)')
ylabel('$(\frac{kgm^2}{s^2})$','Interpreter','latex')
grid on

subplot(2,2,2)
plot(tspan2,u_record2(:,1),'b','LineWidth',2);
hold on
plot(tspan2,u_record2(:,2),'c','LineWidth',2);
plot(tspan2,u_record2(:,3),'m','LineWidth',2);
plot(tspan2,utrans_min2*ones(1,length(tspan2)),'r--','LineWidth',2)
plot(tspan2,utrans_max2*ones(1,length(tspan2)),'r--','LineWidth',2)
hold off
ylim([utrans_min2 -utrans_min2])
legend('$u_x$','$u_y$','$u_z$','constraint','Interpreter','latex')
ylabel('(km/s)','Interpreter','latex')
title('Leg 2')
grid on

subplot(2,2,4)
plot(tspan2,u_record2(:,4),'b','LineWidth',2);
hold on
plot(tspan2,u_record2(:,5),'c','LineWidth',2);
plot(tspan2,u_record2(:,6),'m','LineWidth',2);
plot(tspan2,urot_min2*ones(1,length(tspan2)),'r--','LineWidth',2)
plot(tspan2,urot_max2*ones(1,length(tspan2)),'r--','LineWidth',2)
hold off
ylim([urot_min2 -urot_min2])
legend('$M_1$','$M_2$','$M_3$','constraint','Interpreter','latex')
xlabel('time (s)')
ylabel('$(\frac{kgm^2}{s^2})$','Interpreter','latex')
grid on
%%

% Plot state errors

figure
subplot(4,3,1);
plot(tspan,abs(xi_record1(:,1)-xi_nonvis_record1(:,1)),'b','LineWidth',2);
hold on
plot(tspan(end)+tspan,abs(xi_record2(:,1)-xi_nonvis_record2(:,1)),'b','LineWidth',2);
plot([tspan(end) tspan(end)],...
    [0 max(max(abs(xi_record1(:,1)-xi_nonvis_record1(:,1))),...
    max(abs(xi_record2(:,1)-xi_nonvis_record2(:,1))))],'g','LineWidth',3)
hold off
ylabel('$x$ (km)','Interpreter','latex')
sgtitle('Estimate Error','FontSize',20)
xlim([0 80])
grid on

subplot(4,3,2);
plot(tspan,abs(xi_record1(:,2)-xi_nonvis_record1(:,2)),'b','LineWidth',2);
hold on
plot(tspan(end)+tspan,abs(xi_record2(:,2)-xi_nonvis_record2(:,2)),'b','LineWidth',2);
plot([tspan(end) tspan(end)],...
    [0 max(max(abs(xi_record1(:,2)-xi_nonvis_record1(:,2))),...
    max(abs(xi_record2(:,2)-xi_nonvis_record2(:,2))))],'g','LineWidth',3)
hold off
ylabel('$y$ (km)','Interpreter','latex')
xlim([0 80])
grid on

subplot(4,3,3);
plot(tspan,abs(xi_record1(:,3)-xi_nonvis_record1(:,3)),'b','LineWidth',2);
hold on
plot(tspan(end)+tspan,abs(xi_record2(:,3)-xi_nonvis_record2(:,3)),'b','LineWidth',2);
plot([tspan(end) tspan(end)],...
    [0 max(max(abs(xi_record1(:,3)-xi_nonvis_record1(:,3))),...
    max(abs(xi_record2(:,3)-xi_nonvis_record2(:,3))))],'g','LineWidth',3)
hold off
ylabel('$z$ (km)','Interpreter','latex')
xlim([0 80])
grid on

subplot(4,3,4);
plot(tspan,abs(xi_record1(:,4)-xi_nonvis_record1(:,4)),'b','LineWidth',2);
hold on
plot(tspan(end)+tspan,abs(xi_record2(:,4)-xi_nonvis_record2(:,4)),'b','LineWidth',2);
plot([tspan(end) tspan(end)],...
    [0 max(max(abs(xi_record1(:,4)-xi_nonvis_record1(:,4))),...
    max(abs(xi_record2(:,4)-xi_nonvis_record2(:,4))))],'g','LineWidth',3)
hold off
ylabel('$\dot{x}$ (km/s)','Interpreter','latex')
xlim([0 80])
grid on

subplot(4,3,5);
plot(tspan,abs(xi_record1(:,5)-xi_nonvis_record1(:,5)),'b','LineWidth',2);
hold on
plot(tspan(end)+tspan,abs(xi_record2(:,5)-xi_nonvis_record2(:,5)),'b','LineWidth',2);
plot([tspan(end) tspan(end)],...
    [0 max(max(abs(xi_record1(:,5)-xi_nonvis_record1(:,5))),...
    max(abs(xi_record2(:,5)-xi_nonvis_record2(:,5))))],'g','LineWidth',3)
hold off
ylabel('$\dot{y}$ (km/s)','Interpreter','latex')
xlim([0 80])
grid on

subplot(4,3,6);
plot(tspan,abs(xi_record1(:,6)-xi_nonvis_record1(:,6)),'b','LineWidth',2);
hold on
plot(tspan(end)+tspan,abs(xi_record2(:,6)-xi_nonvis_record2(:,6)),'b','LineWidth',2);
plot([tspan(end) tspan(end)],...
    [0 max(max(abs(xi_record1(:,6)-xi_nonvis_record1(:,6))),...
    max(abs(xi_record2(:,6)-xi_nonvis_record2(:,6))))],'g','LineWidth',3)
hold off
ylabel('$\dot{z}$ (km/s)','Interpreter','latex')
xlim([0 80])
grid on

subplot(4,3,7);
plot(tspan,abs(xi_record1(:,7)-xi_nonvis_record1(:,7)),'b','LineWidth',2);
hold on
plot(tspan(end)+tspan,abs(xi_record2(:,7)-xi_nonvis_record2(:,7)),'b','LineWidth',2);
plot([tspan(end) tspan(end)],...
    [0 max(max(abs(xi_record1(:,7)-xi_nonvis_record1(:,7))),...
    max(abs(xi_record2(:,7)-xi_nonvis_record2(:,7))))],'g','LineWidth',3)
hold off
ylabel('$\phi$ (rad)','Interpreter','latex')
xlim([0 80])
grid on

subplot(4,3,8);
plot(tspan,abs(xi_record1(:,8)-xi_nonvis_record1(:,8)),'b','LineWidth',2);
hold on
plot(tspan(end)+tspan,abs(xi_record2(:,8)-xi_nonvis_record2(:,8)),'b','LineWidth',2);
plot([tspan(end) tspan(end)],...
    [0 max(max(abs(xi_record1(:,8)-xi_nonvis_record1(:,8))),...
    max(abs(xi_record2(:,8)-xi_nonvis_record2(:,8))))],'g','LineWidth',3)
hold off
ylabel('$\theta$ (rad)','Interpreter','latex')
xlim([0 80])
grid on

subplot(4,3,9);
plot(tspan,abs(xi_record1(:,9)-xi_nonvis_record1(:,9)),'b','LineWidth',2);
hold on
plot(tspan(end)+tspan,abs(xi_record2(:,9)-xi_nonvis_record2(:,9)),'b','LineWidth',2);
plot([tspan(end) tspan(end)],...
    [0 max(max(abs(xi_record1(:,9)-xi_nonvis_record1(:,9))),...
    max(abs(xi_record2(:,9)-xi_nonvis_record2(:,9))))],'g','LineWidth',3)
hold off
ylabel('$\psi$ (rad)','Interpreter','latex')
xlim([0 80])
grid on

subplot(4,3,10);
plot(tspan,abs(xi_record1(:,10)-xi_nonvis_record1(:,10)),'b','LineWidth',2);
hold on
plot(tspan(end)+tspan,abs(xi_record2(:,10)-xi_nonvis_record2(:,10)),'b','LineWidth',2);
plot([tspan(end) tspan(end)],...
    [0 max(max(abs(xi_record1(:,10)-xi_nonvis_record1(:,10))),...
    max(abs(xi_record2(:,10)-xi_nonvis_record2(:,10))))],'g','LineWidth',3)
hold off
xlabel('time (s)')
ylabel('$\omega_x$ (rad/s)','Interpreter','latex')
xlim([0 80])
grid on

subplot(4,3,11);
plot(tspan,abs(xi_record1(:,11)-xi_nonvis_record1(:,11)),'b','LineWidth',2);
hold on
plot(tspan(end)+tspan,abs(xi_record2(:,11)-xi_nonvis_record2(:,11)),'b','LineWidth',2);
plot([tspan(end) tspan(end)],...
    [0 max(max(abs(xi_record1(:,11)-xi_nonvis_record1(:,11))),...
    max(abs(xi_record2(:,11)-xi_nonvis_record2(:,11))))],'g','LineWidth',3)
hold off
xlabel('time (s)')
ylabel('$\omega_y$ (rad/s)','Interpreter','latex')
xlim([0 80])
grid on

subplot(4,3,12);
plot(tspan,abs(xi_record1(:,12)-xi_nonvis_record1(:,12)),'b','LineWidth',2);
hold on
plot(tspan(end)+tspan,abs(xi_record2(:,12)-xi_nonvis_record2(:,12)),'b','LineWidth',2);
plot([tspan(end) tspan(end)],...
    [0 max(max(abs(xi_record1(:,12)-xi_nonvis_record1(:,12))),...
    max(abs(xi_record2(:,12)-xi_nonvis_record2(:,12))))],'g','LineWidth',3)
hold off
xlabel('time (s)')
ylabel('$\omega_z$ (rad/s)','Interpreter','latex')
xlim([0 80])
grid on
%%
% Plot histories y
figure
subplot(1,2,1)
plot(y_record1(:,5),y_record1(:,6),'r','LineWidth',2,'HandleVisibility','off')
hold on
plot(y_record1(:,7),y_record1(:,8),'r','LineWidth',2,'HandleVisibility','off')
plot(y_record1(:,9),y_record1(:,10),'r','LineWidth',2,'HandleVisibility','off')
plot(y_record1(:,11),y_record1(:,12),'r','LineWidth',2,'HandleVisibility','off')
plot(y_record1(1,5),y_record1(1,6),'.r','MarkerSize',50,'LineWidth',2)
plot(y_record1(1,7),y_record1(1,8),'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot(y_record1(1,9),y_record1(1,10),'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot(y_record1(1,11),y_record1(1,12),'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot([radius_cnstr1 -radius_cnstr1 -radius_cnstr1 radius_cnstr1 radius_cnstr1],...
    [radius_cnstr1 radius_cnstr1 -radius_cnstr1 -radius_cnstr1 radius_cnstr1],...
    'r--','LineWidth',2);
hold off
xlabel('pixel x')
ylabel('pixel y')
legend('Starting Point',...
    'Constraints','location','Southeast')
axis equal
title('Leg 1')
grid on

subplot(1,2,2)
plot(y_record2(:,5),y_record2(:,6),'r','LineWidth',2,'HandleVisibility','off')
hold on
plot(y_record2(:,7),y_record2(:,8),'r','LineWidth',2,'HandleVisibility','off')
plot(y_record2(:,9),y_record2(:,10),'r','LineWidth',2,'HandleVisibility','off')
plot(y_record2(:,11),y_record2(:,12),'r','LineWidth',2,'HandleVisibility','off')
plot(y_record2(1,5),y_record2(1,6),'.r','MarkerSize',50,'LineWidth',2)
plot(y_record2(1,7),y_record2(1,8),'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot(y_record2(1,9),y_record2(1,10),'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot(y_record2(1,11),y_record2(1,12),'.r','MarkerSize',50,'LineWidth',2,'HandleVisibility','off')
plot([radius_cnstr2 -radius_cnstr2 -radius_cnstr2 radius_cnstr2 radius_cnstr2],...
    [radius_cnstr2 radius_cnstr2 -radius_cnstr2 -radius_cnstr2 radius_cnstr2],...
    'r--','LineWidth',2);
hold off
xlabel('pixel x')
ylabel('pixel y')
legend('Starting Point','location','Southeast')
axis equal
title('Leg 2')
sgtitle('History of $\mathbf{c_1}...\mathbf{c_4}$','Interpreter','latex','FontSize',20)
grid on







