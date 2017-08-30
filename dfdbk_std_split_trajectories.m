% Plot three trajectories of the derivative feedback network with STD on
% the excitatory to excitatory projections.

w = 100;
k = 1.1;
q = 0.5;
diff = -0.005;
tr_e = 0.5;
U_e = 0.1;
delta_U = [0.0,0.0,0.0];
delta_tr = [0.0,0.1,0.25];

delta_tr_leg = cell(length(delta_tr),1);
for i = 1:length(delta_tr)
   delta_tr_leg{i} = ['\Delta \tau_r=',num2str(delta_tr(i))];
end


I_step = 10;
w_in = 1;
q_in = 0.0;
Io_end = 2.0;
t_end = 4.0;

solutions = cell(length(delta_U),1);
eigenvalues = cell(length(delta_U),1);
nonlin_steady_states = cell(length(delta_U),1);

lin_sim = cell(length(delta_U),1);

lin_steady_states = cell(length(delta_U),1);
lin_approx_rise_times = cell(length(delta_U),1);
lin_rise_times = cell(length(delta_U),1);

f = @(x1,x2) I_step./(1-(1-q-diff-(k*w)./(1+k*w)*q).*w.*x1-(q+diff-(k*w)./(1+k*w)*(1-q)).*w.*x2);
rise_times = zeros(length(delta_U),1);
decay_times = zeros(length(delta_U),1);

for i = 1:length(delta_U)
    tr1_e = tr_e + delta_tr(i);
    tr2_e = tr_e - delta_tr(i);
    
    U1_e = U_e + delta_U(i);
    U2_e = U_e - delta_U(i);
    
    [sol,max_eigs,ss_all] = dfdbk_std_split_sim(w,k,q,diff,tr1_e,tr2_e,U1_e,U2_e,I_step,w_in,q_in,Io_end,t_end);
    [rise_time,decay_time,nix] = dfdbk_std_split_RandD(w,k,q,diff,tr1_e,tr2_e,U1_e,U2_e,I_step,w_in,q_in);
    
    solutions{i} = sol;
    eigenvalues{i} = max_eigs;
    
    -log(9)./max_eigs;
    
    nonlin_steady_states{i} = ss_all;
    rise_times(i) = rise_time;
    decay_times(i) = decay_time;
    
    x1  = sol.y(11,:);
    x2  = sol.y(12,:);
    lin_steady_states{i} = f(x1,x2);
    
    
    RT = zeros(length(x1),1);
    for j = 1:length(x1);
        RT(j) = dfdbk_compute_RT_split(x1(j),x2(j),w,q,q_in,diff,k);
    end
    lin_rise_times{i} = RT;
    
    lin_sim{i} = dfdbk_compute_sim(w,k,q,diff,I_step,w_in,q_in,Io_end,t_end);
    
end

% figure('color','w')
% for i = 1:length(delta_U)
%     holder = solutions{i};
%     ss = nonlin_steady_states{i};
%     T = holder.x;
%     ye = holder.y(1,:);
%     plot(T,ye,'LineWidth',2)
%     hold on
%     plot([T(1),T(end)],[ss(1),ss(1)],'k--')
% end
% set(gca,'FontSize',14)
% xlabel('Time (s)','FontSize',25)
% ylabel('Activity (Hz)','FontSize',25)

figure('color','w')
subplot(2,2,1)
for i = 1:length(delta_U)
    holder = lin_sim{i};
    T = holder.x;
    ye = holder.y(1,:);
    plot(T,ye,'LineWidth',2)
    hold on
end
set(gca,'FontSize',14)
xlabel('Time (s)','FontSize',25)
ylabel('Activity (Hz)','FontSize',25)


subplot(2,2,2)
for i = 1:length(delta_U)
    holder = solutions{i};
    T = holder.x;
    ye = holder.y(1,:);
    plot(T,ye,'LineWidth',2)
    hold on
end
set(gca,'FontSize',14)
legend(delta_tr_leg{:})
xlabel('Time (s)','FontSize',25)
ylabel('Activity (Hz)','FontSize',25)



subplot(2,2,3)
for i = 1:length(delta_U)
    holder = solutions{i};
    T = holder.x;
    R_ss = lin_steady_states{i};
    plot(T,R_ss./I_step,'LineWidth',2)
    hold on
    
%     holder = nonlin_steady_states{i}./I_step;
%     plot([T(1),T(end)],[holder(1),holder(1)],'k--')
%     hold on
    
end
set(gca,'FontSize',14)
xlabel('Time (s)','FontSize',25)
ylabel('Gain','FontSize',25)


subplot(2,2,4)
for i = 1:length(delta_U)
    holder = solutions{i};
    T = holder.x;
    
    RT = lin_rise_times{i};
    plot(T,RT,'LineWidth',2)
    hold on
    
%     RT = lin_rise_times{i};
%     plot(T,RT)
%     hold on
    
%     holder = eigenvalues{i};
%     holder = -log(9)./holder(1);
%     plot([T(1),T(end)],[holder,holder],'b--')
%     plot([T(1),T(end)],[rise_times(i),rise_times(i)],'k--')
%     plot([T(1),T(end)],[decay_times(i),decay_times(i)],'r--')
%     hold on
    
end
set(gca,'FontSize',14)
xlabel('Time (s)','FontSize',25)
ylabel('Rise Time (s)','FontSize',25)














