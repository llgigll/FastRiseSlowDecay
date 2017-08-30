% This script sets up the plots for figure 1. There are 6 scripts for
% producing the correct matlab files:
% 1) pfdbk_std_impact.m (this script)
% 2) dfdbk_std_impact.m (figure 2)
% 3) dfdbk_std_on_ee_RandD.m (figures 3 and 4)
% 4) explain_split_std.m (figure 5)
% 5) dfdbk_std_manysynapses_sim_ex.m (figures 6 and 8)
% 6) not built yet (figure 7 and spikes plots in figures 1 and 2)

% Show how STD on recurrent connections impacts the response time of the
% positive feedback network.

w = 0.9936;
q = 0.5;
tr_e = 0.5;
U_e_hold = [0.2,0.1,0.05,0.0];
u_leg = {['u=0.2'],['u=0.10'],['u=0.05'],['No STD']};
w_in = 1;
q_in = 0.5;
Io_end = 2.0;
t_end = 4.0;

R_ss = 20;


f = @(xw) I_step./(1-xw);
I = @(R,w,tr,u)(1-w/(1+tr*u*R))*R;
solutions = cell(length(U_e_hold),1);

lin_rise_times = cell(length(U_e_hold),1);


for i = 1:length(U_e_hold)
    U_e = U_e_hold(i);
    if U_e == 0
        I_step = R_ss*(1-w);
        sol = pfdbk_compute_sim(w,q,I_step,w_in,q_in,Io_end,t_end);
        x = ones(length(sol.x),1);
    else
        I_step = I(R_ss,w,tr_e,U_e);
        [sol,ss_all] = pfdbk_std_on_ee_sim(w,q,tr_e,U_e,I_step,w_in,q_in,Io_end,t_end);
        x  = sol.y(6,:);
    end

    solutions{i} = sol;


    RT = zeros(length(x),1);
    for j = 1:length(x);
        if x(j)*w>0.999
            RT(j) = pfdbk_compute_RT(0.99,q,q_in);
        else
            RT(j) = pfdbk_compute_RT(x(j)*w,q,q_in);
        end
    end

    lin_rise_times{i} = RT;
    

end


figure('color','w')
for i = 1:length(U_e_hold)
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
ylim([0,37])
set(gca,'FontSize',20)
xlabel('Time (s)','FontSize',30)
ylabel('Decay time (s)','FontSize',30)
legend(u_leg{:})


U_e_hold = [0.2,0.1,0.05];
w_hold = [0.01:0.01:0.99,1.0];
decay_times = zeros(length(w_hold),length(U_e_hold));
for j = 1:length(U_e_hold)
    U_e = U_e_hold(j);
    for i = 1:length(w_hold)
        w = w_hold(i);
        I_step = I(R_ss,w,tr_e,U_e);    
        [rise_time,decay_time] = pfdbk_std_on_ee_RandD(w,q,tr_e,U_e,I_step,w_in,q_in);
        decay_times(i,j) = decay_time;
    end
end


figure('color','w')
for i = 1:length(U_e_hold)
    plot(w_hold,decay_times(:,i),'LineWidth',2)
    hold on
end
ylim([0,1.25*max(decay_times(:,3))])
set(gca,'FontSize',20)
xlabel('w','FontSize',30)
ylabel('Decay Time (s)','FontSize',30)
legend(u_leg{:})
xlim([0,1])

colors = get(gca,'ColorOrder');
w = 0.01:0.01:0.99;
rise = pfdbk_std_impact_linrisetime(w);
figure('color','w')
plot(w,rise,'color',colors(4,:),'LineWidth',2)
xlabel('w','FontSize',30)
ylabel('Decay time (s)','FontSize',30)
ax = gca; 
ax.FontSize = 20;


U_e_hold = [0.01:0.01:0.5];
w = 1.0;
decay_times = zeros(length(U_e_hold),1);
for j = 1:length(U_e_hold)
    U_e = U_e_hold(j);
    I_step = I(R_ss,w,tr_e,U_e);
    [rise_time,decay_time] = pfdbk_std_on_ee_RandD(w,q,tr_e,U_e,I_step,w_in,q_in);
    decay_times(j) = decay_time;
end


figure('color','w')
plot(U_e_hold,decay_times,'k','LineWidth',2)
set(gca,'FontSize',20)
xlabel('u','FontSize',30)
ylabel('Decay Time (s)','FontSize',30)
xlim([0,0.3])
ylim([0,3.0])





