% This script sets up the plots for figure 2.

% Plot three trajectories of the derivative feedback network with STD on
% the excitatory to excitatory projections.

w = 100;
k = 1.1;
q = 0.5;
diff = 0.115;
U_e_hold = [0.2,0.1,0.05,0.0];
u_leg = {['u=0.20'],['u=0.10'],['u=0.05'],['No STD']};
tr_e = 0.5;
w_in = 1;
q_in = 0.5;

Io_end = 2.0;
t_end = 4.0;
Re_ss = 20;

t_ampa = 0.005;
t_nmda = 0.1;
tau_ie = 0.5*t_ampa+0.5*t_nmda;
tau_ee = (0.5-diff)*t_ampa+(0.5+diff)*t_nmda;
delta_tau = (tau_ee - tau_ie)*10^3;



solutions = cell(length(U_e_hold),1);
lin_rise_times = cell(length(diff),1);

for i = 1:length(U_e_hold)
    U_e = U_e_hold(i);
    if U_e == 0
        I_step = Re_ss*(1-w/(1+k*w));
        sol = dfdbk_compute_sim(w,k,q,diff,I_step,w_in,q_in,Io_end,t_end);
        x = ones(length(sol.x),1);
    else
        I_step = Re_ss*(1-(w/(1+k*w))*(1/(1+tr_e*U_e*Re_ss)));
        [sol,max_eigs,ss_all] = dfdbk_std_on_ee_sim(w,k,q,diff,tr_e,U_e,I_step,w_in,q_in,Io_end,t_end);
        x  = sol.y(11,:);
    end
    
    solutions{i} = sol;

    RT = zeros(length(x),1);
    for j = 1:length(x);
        RT(j) = dfdbk_compute_RT(x(j)*w,q,q_in,diff,k./x(j));
    end
    lin_rise_times{i} = RT;
    
    
end


figure('color','w')
for i = 1:(length(U_e_hold))
    holder = solutions{i};
    T = holder.x;
    ye = holder.y(1,:);
    if U_e_hold(i) == 0.0
        plot(T,ye,'LineWidth',2)
    else
        plot(T,ye,'LineWidth',2)
    end
    hold on
end
set(gca,'FontSize',20)
xlabel('Time (s)','FontSize',30)
ylabel('Activity (Hz)','FontSize',30)
legend(u_leg{:})


figure('color','w')
for i = 1:length(U_e_hold)
    holder = solutions{i};
    T = holder.x;
    if U_e_hold(i) == 0.0
        plot(T,lin_rise_times{i},'LineWidth',2)
    else
        plot(T,lin_rise_times{i},'LineWidth',2)
    end
    hold on
end
set(gca,'FontSize',20)
ylim([0,35])
xlabel('Time (s)','FontSize',30)
ylabel('Decay time (s)','FontSize',30)
legend(u_leg{:})






U_e_hold = [0.2,0.1,0.05];
w_hold = [1:4:200];
decay_times = zeros(length(w_hold),length(U_e_hold));
for j = 1:length(U_e_hold)
    U_e = U_e_hold(j);
    for i = 1:length(w_hold)
        w = w_hold(i);
        I_step = Re_ss*(1-(w/(1+k*w))*(1/(1+tr_e*U_e*Re_ss)));
        [rise_time,decay_time,nix] = dfdbk_std_on_ee_RandD(w,k,q,diff,tr_e,U_e,I_step,w_in,q_in);
        decay_times(i,j) = decay_time;
    end
end


figure('color','w')
for i = 1:length(U_e_hold)
    plot(w_hold,decay_times(:,i),'LineWidth',2)
    hold on
end
ylim([0,1.25*max(decay_times(:,end))])
set(gca,'FontSize',20)
xlabel('w','FontSize',30)
ylabel('Decay Time (s)','FontSize',30)
legend(u_leg{:})
colors = get(gca,'ColorOrder');


w = 1:4:200;
k = 1.1;
delta_q = 0.115;
rise = dfdbk_std_impact_linrisetime(delta_q,w,k);

figure('color','w')
for i = 1:length(k)
    plot(w,rise,'color',colors(4,:),'LineWidth',2)
    hold on
end
ylabel('Decay time (s)','FontSize',30)
xlabel('w','FontSize',30)
ax = gca; 
ax.FontSize = 20;







