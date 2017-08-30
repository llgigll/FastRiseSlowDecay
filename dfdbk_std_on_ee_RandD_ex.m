% This script sets up the plots for figures 3 and 4 showing that changing
% delta tau also significantly increases the rise time and that feedforward
% STD is not sufficient to fix this.


Io_end = 2.0;
t_end = 4.0;
w = 100;
k = 1.1;
q = 0.5;
diff_grid = [0.0:0.0005:0.01]./0.095;
diff_plot = [0,0.005,0.01]./0.095;
U_e = [0.1,0.2,0.3];
tr_e = 0.5;
I_step = 20;
w_in = 0.5;
w_in_base = w_in;
q_in = 0.5;

rise = zeros(length(diff_grid),length(U_e));
decay = zeros(length(diff_grid),length(U_e));
for i = 1:length(diff_grid)
    parfor j = 1:length(U_e)
        [rise_time,decay_time,max_eigs] = dfdbk_std_on_ee_RandD(w,k,q,diff_grid(i),...
                                                        tr_e,U_e(j),I_step,w_in,q_in);
        rise(i,j) = rise_time;
        decay(i,j) = decay_time;
    end
end

t_ampa = 0.005;
t_nmda = 0.1;
tau_ie = 0.5*t_ampa+0.5*t_nmda;
tau_ee = (0.5-diff_grid)*t_ampa+(0.5+diff_grid)*t_nmda;
delta_tau = (tau_ee - tau_ie)*10^3;

U_e_leg = cell(length(U_e),1);
for i = 1:length(U_e)
   U_e_leg{i} = ['u_{e}=',num2str(U_e(i))];
end

figure('color','w')
mylines = {'--k','-.k','-k'};
for l = 1:length(U_e)
    plot(delta_tau,rise(:,l),mylines{l},'LineWidth',2)
    hold on
end
xlabel('\Delta\tau (ms)','FontSize',30)
ylabel('Rise time (s)','FontSize',30)
set(gca,'FontSize',20)
legend(U_e_leg{:},'Location','NorthWest')

figure('color','w')
for l = 1:length(U_e)
    plot(delta_tau,decay(:,l),mylines{l},'LineWidth',2)
    hold on
end
xlabel('\Delta\tau (ms)','FontSize',30)
ylabel('Decay time (s)','FontSize',30)
set(gca,'FontSize',20)
legend(U_e_leg{:},'Location','NorthWest')


figure('color','w')
for l = 1:length(U_e)
    plot(delta_tau,rise(:,l)./decay(:,l),mylines{l},'LineWidth',2)
    hold on
end
ylim([0,0.8])
xlabel('\Delta\tau (ms)','FontSize',30)
ylabel('Rise time / decay time','FontSize',30)
set(gca,'FontSize',20)
legend(U_e_leg{:})



t_ampa = 0.005;
t_nmda = 0.1;
tau_ie = 0.5*t_ampa+0.5*t_nmda;
tau_ee = (0.5-diff_plot)*t_ampa+(0.5+diff_plot)*t_nmda;
delta_tau = (tau_ee - tau_ie)*10^3;

delta_tau_leg = cell(length(delta_tau),1);
for i = 1:length(delta_tau)
   delta_tau_leg{i} = ['\Delta\tau=',num2str(delta_tau(i)),' ms'];
end

solutions = cell(length(diff_plot),1);

for i = 1:length(diff_plot)
    [sol,max_eigs,ss_all] = dfdbk_std_on_ee_sim(w,k,q,diff_plot(i),tr_e,U_e(1),I_step,w_in,q_in,Io_end,t_end);    
    solutions{i} = sol;
    
end

figure('color','w')
for i = 1:length(diff_plot)
    holder = solutions{i};
    T = holder.x;
    ye = holder.y(1,:);
    plot(T,ye,'LineWidth',2)
    hold on
end
set(gca,'FontSize',20)
legend(delta_tau_leg{:})
xlabel('Time (s)','FontSize',30)
ylabel('Activity (Hz)','FontSize',30)





 


q_in = 0.5;
diff = interp1(decay(:,1),diff_grid,1.0);
% t_ampa = 0.005;
% t_nmda = 0.1;
% tau_ie = 0.5*t_ampa+0.5*t_nmda;
% tau_ee = (0.5-diff)*t_ampa+(0.5+diff)*t_nmda;
% delta_tau = (tau_ee - tau_ie)*10^3;
U_e = 0.1;
tr_ff = 0.5;
U_ff = fliplr([0.1,0.25,0.5]);
solutions = cell(length(U_ff),1);
w_in_hold = zeros(length(U_ff),1);
for i = 1:length(U_ff)
    w_in = w_in_base*(1+tr_ff*U_ff(i)*I_step);
    w_in_hold(i) = w_in;
    [sol,ss_all] = dfdbk_std_on_ee_eo_sim(w,k,q,diff,tr_e,U_e,tr_ff,U_ff(i),I_step,w_in,q_in,Io_end,t_end);
    solutions{i} = sol;
end

U_ff_leg = cell(length(U_ff),1);
for i = 1:length(delta_tau)
   U_ff_leg{i} = ['u_{ff}=',num2str(U_ff(i))];
end

figure('color','w')
for i = 1:length(U_ff)
    T = solutions{i}.x;
    ye = solutions{i}.y(1,:);
    plot(T,ye,'LineWidth',2)
    hold on
end
xlabel('Time (s)','FontSize',30)
ylabel('Activity (Hz)','FontSize',30)
set(gca,'FontSize',20)
legend(U_ff_leg{:})

figure('color','w')
for i = 1:length(U_ff)
    T = solutions{i}.x;
    x_ff = solutions{i}.y(end,:);
    plot(T,x_ff,'LineWidth',2)
    hold on
end
xlabel('Time (s)','FontSize',30)
ylabel('STD','FontSize',30)
set(gca,'FontSize',20)

figure('color','w')
for i = 1:length(U_ff)
    T = solutions{i}.x;
    Sa_in = solutions{i}.y(9,:);
    Sn_in = solutions{i}.y(10,:);
    I_in = w_in_hold(i)*((1-q_in)*Sa_in+q_in*Sn_in);
    plot(T,I_in,'LineWidth',2)
    hold on
end
xlabel('Time (s)','FontSize',30)
ylabel('I_{in}','FontSize',30)
set(gca,'FontSize',20)
legend(U_ff_leg{:})




t_input = -0.5:0.001:t_end;
input = (t_input>0).*(t_input<Io_end)*I_step;
figure('color','w')
plot(t_input,input,'color',[0 .5 0],'LineWidth',2)
ylim([0,I_step*1.2])
xlabel('Time (s)','FontSize',30)
ylabel('R_{in} (Hz)','FontSize',30)
set(gca,'FontSize',20)












