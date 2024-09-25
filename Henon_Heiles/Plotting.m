% Plotting

% all
% figure(1)
% subplot(1,3,1); plot(omega,abs(FT3_ra))
% legend('$S(\omega )$','Interpreter','latex','FontSize',20,'Linewidth',2);
% title('Rank-adaptive BUG')
% subplot(1,3,2); plot(omega,abs(FT3_BUG))
% legend('$S(\omega )$','Interpreter','latex','FontSize',20,'Linewidth',2)
% title('Fixed-rank BUG')
% subplot(1,3,3); plot(omega,abs(FT3))
% legend('$S(\omega )$','Interpreter','latex','FontSize',20,'Linewidth',2)
% title('Parallel BUG')
% 
% figure(2)
% subplot(2,3,1); plot(0:dt:T_end,abs(Norm_con_ad))
% title('Rank-adaptive BUG')
% legend('$ \vert \vert X(t) \vert \vert$','Interpreter','latex','FontSize',14)
% subplot(2,3,4); plot(0:dt:T_end,abs(Energy_con_ad))
% legend('$ < X(t) \vert H \vert X(t) > $','Interpreter','latex','FontSize',14)
% 
% subplot(2,3,2); plot(0:dt:T_end,abs(Norm_con_BUG))
% title('Fixed-rank BUG')
% legend('$ \vert \vert X(t) \vert \vert$','Interpreter','latex','FontSize',14)
% subplot(2,3,5); plot(0:dt:T_end,abs(Energy_con_BUG))
% legend('$ < X(t) \vert H \vert X(t) > $','Interpreter','latex','FontSize',14)
% 
% subplot(2,3,3); plot(0:dt:T_end,abs(Norm_con))
% title('Parallel BUG')
% legend('$ \vert \vert X(t) \vert \vert$','Interpreter','latex','FontSize',14)
% subplot(2,3,6); plot(0:dt:T_end,abs(Energy_con))
% legend('$ < X(t) \vert H \vert X(t) > $','Interpreter','latex','FontSize',14)
% 
% figure(3)
% plot(0:dt:T_end,rank_max_ad,'Linewidth',2)
% hold on
% plot(0:dt:T_end,rank_max_BUG,'Linewidth',2)
% plot(0:dt:T_end,rank_max,'Linewidth',2)
% legend('Rank-adaptive BUG','Fixed-rank BUG','Parallel BUG','FontSize', 20)

% without parallel plus autocorrelation function
load('HH_d8_Tend30_dt0005_rmax4.mat')
figure(1)
subplot(1,2,1); plot(omega,abs(FT3_ra),'Linewidth',1)
legend('$S(\omega )$','Interpreter','latex','FontSize',20,'Linewidth',2);
title('Rank-adaptive BUG')
subplot(1,2,2); plot(omega,abs(FT3_BUG),'Linewidth',1)
legend('$S(\omega )$','Interpreter','latex','FontSize',20,'Linewidth',2)
title('Fixed-rank BUG')

figure(2)
subplot(2,2,1); plot(0:dt:T_end,abs(Norm_con_ad),'Linewidth',2)
title('Rank-adaptive BUG')
legend('$ \vert \vert X(t) \vert \vert$','Interpreter','latex','FontSize',14)
subplot(2,2,3); plot(0:dt:T_end,abs(Energy_con_ad),'Linewidth',2)
legend('$ < X(t) \vert H \vert X(t) > $','Interpreter','latex','FontSize',14)

subplot(2,2,2); plot(0:dt:T_end,abs(Norm_con_BUG),'Linewidth',2)
title('Fixed-rank BUG')
legend('$ \vert \vert X(t) \vert \vert$','Interpreter','latex','FontSize',14)
subplot(2,2,4); plot(0:dt:T_end,abs(Energy_con_BUG),'Linewidth',2)
legend('$ < X(t) \vert H \vert X(t) > $','Interpreter','latex','FontSize',14)

figure(3)
subplot(2,1,1); plot(0:dt:T_end,ac_ad,'Linewidth',2)
legend('a(t)')
title('Rank-adaptive BUG')
hold on
subplot(2,1,2); plot(0:dt:T_end,ac_BUG,'Linewidth',2)
legend('a(t)')
title('Fixed-rank BUG')

% figure(3)
% plot(0:dt:T_end,rank_max_ad,'Linewidth',2)
% hold on
% plot(0:dt:T_end,rank_max_BUG,'Linewidth',2)
% legend('Rank-adaptive BUG','Fixed-rank BUG','FontSize', 20)

%% plot for results of only fixed-rank
clear all;
load('HH_d8_Tend30_dt0005_rmax4.mat')
subplot(1,3,1); plot(omega,abs(FT3_BUG),'Linewidth',1)
legend('$S(\omega )$','Interpreter','latex','FontSize',20,'Linewidth',2)
title('8 particles')
clear all;
load('HH_d10_Tend60_dt0005_rmax4.mat')
subplot(1,3,2); plot(omega,abs(FT3_BUG),'Linewidth',1)
legend('$S(\omega )$','Interpreter','latex','FontSize',20,'Linewidth',2)
title('10 particles')
clear all;
load('HH_d16_Tend50_dt0005_rmax6.mat')
subplot(1,3,3); plot(omega,abs(FT3_BUG),'Linewidth',1)
legend('$S(\omega )$','Interpreter','latex','FontSize',20,'Linewidth',2)
title('16 particles')

