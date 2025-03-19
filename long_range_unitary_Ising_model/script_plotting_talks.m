% Plotting for Talks

figure(1)
subplot(1,2,1)
plot([0 time],abs(nn_ad),'*-','Linewidth',2)
hold on
plot([0 time],abs(nn_BUG),'x-','Linewidth',2)
plot([0 time],abs(nn_par),'o-','Linewidth',2)
xlabel('Time')
title('Norm')
legend('Rank-adaptive BUG','Fixed-rank BUG','Parallel BUG')

subplot(1,2,2)
plot([0 time],abs(en_ad),'*-','Linewidth',2)
hold on
plot([0 time],abs(en_BUG),'x-','Linewidth',2)
plot([0 time],abs(en_par),'o-','Linewidth',2)
xlabel('Time')
title('Energy')
legend('Rank-adaptive BUG','Fixed-rank BUG','Parallel BUG')


figure(10)
subplot(1,2,1)
plot([0 time],abs(nn_ad),'*-','Linewidth',2)
hold on
plot([0 time],abs(nn_BUG),'x-','Linewidth',2)
plot([0 time],abs(nn_par),'o-','Linewidth',2)
xlabel('Time')
title('Norm')
legend('Rank-adaptive BUG','Fixed-rank BUG','Parallel BUG')

subplot(1,2,2)
semilogy([0 time],abs(abs(nn_ad)- nn_ad(1)),'*-','Linewidth',2)
hold on
semilogy([0 time],abs(abs(nn_BUG)- nn_ad(1)),'x-','Linewidth',2)
semilogy([0 time],abs(abs(nn_par)- nn_ad(1)),'o-','Linewidth',2)
xlabel('Time')
title('Error Norm')
legend('Rank-adaptive BUG','Fixed-rank BUG','Parallel BUG')

figure(11)
subplot(1,2,1)
plot([0 time],abs(en_ad),'*-','Linewidth',2)
hold on
plot([0 time],abs(en_BUG),'x-','Linewidth',2)
plot([0 time],abs(en_par),'o-','Linewidth',2)
xlabel('Time')
title('Energy')
legend('Rank-adaptive BUG','Fixed-rank BUG','Parallel BUG')

subplot(1,2,2)
semilogy(time,abs(abs(en_ad(2:end))- en_ad(1)),'*-','Linewidth',2)
hold on
semilogy(time,abs(abs(en_BUG(2:end))- en_ad(1)),'x-','Linewidth',2)
semilogy(time,abs(abs(en_par(2:end))- en_ad(1)),'o-','Linewidth',2)
xlabel('Time')
title('Error Energy')
legend('Rank-adaptive BUG','Fixed-rank BUG','Parallel BUG')