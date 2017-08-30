% Script for plots showing how distributions of STD can change the average
% time constant of a projection.  This is figure 5.

std = @(t,tr,u,R)1./(1+tr.*u.*R)+(1-(1./(1+tr.*u.*R))).*exp(-(1./tr+u.*R).*t);

u = 0.15;
R = 20;
tr = 0.5;
t = 0:0.01:1.0;
percent = 0.15;

x2 = std(t,tr*(1+percent),u*(1+percent),R);
x1 = std(t,tr*(1-percent),u*(1-percent),R);

figure('color','w')
plot(t,x1,'LineWidth',2)
hold on
plot(t,x2,'LineWidth',2)
set(gca,'FontSize',20)
ylim([0,1])
legend('Weak STD','Strong STD')
xlabel('Time (s)','FontSize',30)
ylabel('STD','FontSize',30)

figure('color','w')
plot(t,(x1+x2)/2,'LineWidth',2)
set(gca,'FontSize',20)
ylim([0,1])
xlabel('Time (s)','FontSize',30)
ylabel('Synaptic Strength','FontSize',30)

tau1 = 0.25*0.005+0.75*0.1;
tau2 = 0.75*0.005+0.25*0.1;
figure('color','w')
plot(t,1000*(x1*tau1+x2*tau2)./(x1+x2),'LineWidth',2)
hold on
plot(t,1000*(x2*tau1+x1*tau2)./(x1+x2),'LineWidth',2)
set(gca,'FontSize',20)
legend('EE','IE')
xlabel('Time (s)','FontSize',30)
ylabel('Average \tau (ms)','FontSize',30)
