

figure; 
x = linspace(0,numel(nDataMean)/480,numel(nDataMean));
subplot(4,1,1); plot(x,nDataMean/1000)
ylabel('Position (m)'); 

x = linspace(0,numel(nVel)/480,numel(nVel));
subplot(4,1,2); plot(x,nVel/1000)
ylabel('Velocity (m/s)');

x = linspace(0,numel(nJerk)/480,numel(nAccel));
subplot(4,1,3); plot(x,nAccel/1000)
hold on; xline((locs(1)+200)/480); xline((locs(2)+200)/480)
ylabel('Acc (m/s^2)'); xlabel('time (s)'); 

x = linspace(0,numel(nJerk)/480,numel(nJerk));
subplot(4,1,4); plot(x,nJerk/1000)
hold on; xline((locs(1)+200)/480); xline((locs(2)+200)/480)
ylabel('Jerk (m/s^3)'); xlabel('time (s)'); 
