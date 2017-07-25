%load plottingdata.dat;
figure(1);
plot(time, x, 'r');
title('Particle Filter Estimation');
xlabel('timesteps'), ylabel('displacement');
hold on;
plot(time, x_hat, 'g');
legend('Truth x','Estimate x','Location','NorthEast')

figure(2);
plot(time, dx, 'r');
hold on;
plot(time, dx_hat, 'g');
title('Particle Filter Estimation');
xlabel('timesteps'), ylabel('velocity');
legend('Truth dx','Estimate dx','Location','NorthEast');