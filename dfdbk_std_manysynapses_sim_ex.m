% Produces plots for figures 6 and 8.

% Plot three trajectories of the derivative feedback network with STD on
% the excitatory to excitatory projections.

I_step = 15.0;

w = 100;
k = 1.1;

w_in = 1.0;
q_in = 0.5;


q = [0.25,0.75]';
U_base = 0.1;
tr_base = 0.5;
qdiff = -0.0075;
diff_template = [ones(1,1);-ones(1,1)];

stim = 2.0;
isi = 2.0;
num_reps = 1;
stim_setup = [stim,isi,num_reps];

%%%%  Simulations for three different percent values.
percent_hold = [0.05,0.1,0.15];
solutions = cell(length(percent_hold),1);
for i = 1:length(percent_hold)
    percent = percent_hold(i);
    delta_U = percent*U_base;
    delta_tr = percent*tr_base;

    U = U_base*ones(size(q));
    U = U+delta_U*diff_template;
    
    tr = tr_base*ones(size(q));
    tr = tr+delta_tr*diff_template;
    
    [sol, ss_all] = dfdbk_std_manysynapses_sim(w,k,q,qdiff,tr,U,I_step,w_in,q_in,stim_setup);
    solutions{i} = sol;
end

percent_leg = cell(length(percent_hold),1);
for i = 1:length(percent_hold)
   percent_leg{i} = ['p=',num2str(percent_hold(i))];
end

figure('color','w')
for i = 1:length(percent_hold)
    holder = solutions{i};
    T = holder.x;
    ye = holder.y(1,:);
    plot(T,ye,'LineWidth',2)
    hold on
end
set(gca,'FontSize',20)
legend(percent_leg{:})
xlabel('Time (s)','FontSize',30)
ylabel('Activity (Hz)','FontSize',30)


figure('color','w')
for i = 1:length(percent_hold)
    holder = solutions{i};
    x1 = holder.y(end-1,:);
    x2 = holder.y(end,:);
    delta_tau = (q(2)-q(1))*(0.005-0.1)*(x1-x2)./(x1+x2)+qdiff*(0.1-0.005);
    T = holder.x;
    plot(T,delta_tau*1000,'LineWidth',2)
    hold on
end
set(gca,'FontSize',20)
legend(percent_leg{:})
xlabel('Time (s)','FontSize',30)
ylabel('\Delta\tau (ms)','FontSize',30)

%%%% Rise and decay across many percent values.

% percent_hold = [0:0.005:0.15];
% rise_all = zeros(length(percent_hold),1);
% decay_all = zeros(length(percent_hold),1);
% 
% parfor i = 1:length(percent_hold)
%     percent = percent_hold(i);
%     delta_U = percent*U_base;
%     delta_tr = percent*tr_base;
% 
%     U = U_base*ones(size(q));
%     U = U+delta_U*diff_template;
%     
%     tr = tr_base*ones(size(q));
%     tr = tr+delta_tr*diff_template;
%     
%     [rise_time, decay_time] = dfdbk_std_manysynapses_RandD(w,k,q,qdiff,tr,U,I_step,w_in,q_in,stim_setup);
%     rise_all(i) = rise_time;
%     decay_all(i) = decay_time;
% end
% 
% figure('color','w')
% plot(percent_hold,rise_all,'--k','LineWidth',2)
% hold on
% plot(percent_hold,decay_all,'k','LineWidth',2)
% 
% set(gca,'FontSize',20)
% legend('Rise time','Decay time')
% xlabel('p','FontSize',30)
% ylabel('Response Time (s)','FontSize',30)
% 
% 
% 
% 
% %%%% Rise and decay across many base values of u
% 
% percent = 0.1;
% U_base_hold = [0.05:0.005:0.4];
% rise_all = zeros(length(percent_hold),1);
% decay_all = zeros(length(percent_hold),1);
% 
% parfor i = 1:length(U_base_hold)
%     U_base = U_base_hold(i);
%     delta_U = percent*U_base;
%     delta_tr = percent*tr_base;
% 
%     U = U_base*ones(size(q));
%     U = U+delta_U*diff_template;
%     
%     tr = tr_base*ones(size(q));
%     tr = tr+delta_tr*diff_template;
%     
%     [rise_time, decay_time] = dfdbk_std_manysynapses_RandD(w,k,q,qdiff,tr,U,I_step,w_in,q_in,stim_setup);
%     rise_all(i) = rise_time;
%     decay_all(i) = decay_time;
% end
% 
% figure('color','w')
% plot(U_base_hold,rise_all,'--k','LineWidth',2)
% hold on
% plot(U_base_hold,decay_all,'k','LineWidth',2)
% 
% set(gca,'FontSize',20)
% legend('Rise time','Decay time')
% xlabel('u','FontSize',30)
% ylabel('Response Time (s)','FontSize',30)
% 
% 
% 
% 
% %%%% Simulations for three different percent values but with FFD STD
% tr_ff = 0.5;
% U_ff = 0.2;
% U_base = 0.1;
% qdiff = 0.0;
% w_in = 1+tr_ff*U_ff*I_step;
% 
% percent_hold = [0.0,0.05,0.1,0.15];
% solutions = cell(length(percent_hold),1);
% for i = 1:length(percent_hold)
%     percent = percent_hold(i);
%     delta_U = percent*U_base;
%     delta_tr = percent*tr_base;
% 
%     U = U_base*ones(size(q));
%     U = U+delta_U*diff_template;
%     
%     tr = tr_base*ones(size(q));
%     tr = tr+delta_tr*diff_template;
%     [sol, ss_all] = dfdbk_ffstd_manysynapses_sim(w,k,q,qdiff,tr,U,tr_ff,U_ff,I_step,w_in,q_in,stim_setup);
%     solutions{i} = sol;
% end
% 
% percent_leg = cell(length(percent_hold),1);
% for i = 1:length(percent_hold)
%    percent_leg{i} = ['p=',num2str(percent_hold(i))];
% end
% 
% figure('color','w')
% for i = 1:length(percent_hold)
%     holder = solutions{i};
%     T = holder.x;
%     ye = holder.y(1,:);
%     plot(T,ye,'LineWidth',2)
%     hold on
% end
% set(gca,'FontSize',20)
% legend(percent_leg{:})
% xlabel('Time (s)','FontSize',30)
% ylabel('Activity (Hz)','FontSize',30)
% 
% 
% 
% %%%% Rise and decay across many percent values but with FFD STD
% 
% 
% percent_hold = [0:0.005:0.15];
% rise_all = zeros(length(percent_hold),1);
% decay_all = zeros(length(percent_hold),1);
% 
% parfor i = 1:length(percent_hold)
%     percent = percent_hold(i);
%     delta_U = percent*U_base;
%     delta_tr = percent*tr_base;
% 
%     U = U_base*ones(size(q));
%     U = U+delta_U*diff_template;
%     
%     tr = tr_base*ones(size(q));
%     tr = tr+delta_tr*diff_template;
%     
%     [rise_time, decay_time] = dfdbk_ffstd_manysynapses_RandD(w,k,q,qdiff,tr,U,tr_ff,U_ff,I_step,w_in,q_in,stim_setup);
%     rise_all(i) = rise_time;
%     decay_all(i) = decay_time;
% end
% 
% figure('color','w')
% plot(percent_hold,rise_all,'--k','LineWidth',2)
% hold on
% plot(percent_hold,decay_all,'k','LineWidth',2)
% 
% set(gca,'FontSize',20)
% legend('Rise time','Decay time')
% xlabel('p','FontSize',30)
% ylabel('Response Time (s)','FontSize',30)
% 
% 
% 
% %%%% Plot three simulations as a function of presentation length
% 
% 
% tr_ff = 0.5;
% U_ff = 0.2;
% U_base = 0.1;
% qdiff = 0.0;
% w_in = 1+tr_ff*U_ff*I_step;
% 
% percent = 0.1;
% delta_U = percent*U_base;
% delta_tr = percent*tr_base;
% 
% U = U_base*ones(size(q));
% U = U+delta_U*diff_template;
% 
% tr = tr_base*ones(size(q));
% tr = tr+delta_tr*diff_template;
% 
% 
% stim_hold = [0.1,0.25,0.5];
% isi = 2.0;
% num_reps = 1;
% 
% solutions = cell(length(stim_hold),1);
% for i = 1:length(stim_hold)
%     stim_setup = [stim_hold(i),isi,num_reps];
%     [sol, ss_all] = dfdbk_ffstd_manysynapses_sim(w,k,q,qdiff,tr,U,tr_ff,U_ff,I_step,w_in,q_in,stim_setup);
%     solutions{i} = sol;
% end
% 
% percent_leg = cell(length(stim_hold),1);
% for i = 1:length(stim_hold)
%    stim_leg{i} = ['Stimulus Length=',num2str(stim_hold(i)),' s'];
% end
% 
% figure('color','w')
% for i = 1:length(stim_hold)
%     holder = solutions{i};
%     T = holder.x;
%     ye = holder.y(1,:);
%     plot(T-stim_hold(i),ye,'LineWidth',2)
%     hold on
% end
% set(gca,'FontSize',20)
% legend(stim_leg{:})
% xlabel('Time (s)','FontSize',30)
% ylabel('Activity (Hz)','FontSize',30)
% 
% 
% 
% %%%% Plot decay time as a function of presentation length
% 
% stim_hold = [0.1:0.05:1.0];
% isi = 4.0;
% num_reps = 1;
% 
% 
% decay_all = zeros(length(stim_hold),1);
% parfor i = 1:length(stim_hold)
%     stim_setup = [stim_hold(i),isi,num_reps];
%     [sol, ss_all] = dfdbk_ffstd_manysynapses_sim(w,k,q,qdiff,tr,U,tr_ff,U_ff,I_step,w_in,q_in,stim_setup);
%     
%     T = sol.x;
%     [na,ind] = find(T<stim_hold(i),1,'Last');
%     ye = sol.y(1,:);
%     ss_val = ye(ind);
%     [na,start_ind] = find(ye>(0.9*ss_val),1,'Last');
%     [na,end_ind] = find(ye>(0.1*ss_val),1,'Last');
%     decay_all(i) = T(end_ind)-T(start_ind);
% end
% 
% 
% tr_ff = 0.5;
% U_ff = 0.2;
% U_base = 0.1;
% w_in = 1+tr_ff*U_ff*I_step;
% qdiff = 0.05;
% 
% percent = 0.0;
% delta_U = percent*U_base;
% delta_tr = percent*tr_base;
% 
% U_base = 0.1;
% U = U_base*ones(size(q));
% U = U+delta_U*diff_template;
% 
% tr_base = 0.5;
% tr = tr_base*ones(size(q));
% tr = tr+delta_tr*diff_template;
% 
% 
% decay_qd_all = zeros(length(stim_hold),1);
% 
% parfor i = 1:length(stim_hold)
%     stim_setup = [stim_hold(i),isi,num_reps];
%     [sol, ss_all] = dfdbk_ffstd_manysynapses_sim(w,k,q,qdiff,tr,U,tr_ff,U_ff,I_step,w_in,q_in,stim_setup);
%     T = sol.x;
%     [na,ind] = find(T<stim_hold(i),1,'Last');
%     ye = sol.y(1,:);
%     ss_val = ye(ind);
%     [na,start_ind] = find(ye>(0.9*ss_val),1,'Last');
%     [na,end_ind] = find(ye>(0.1*ss_val),1,'Last');
%     decay_qd_all(i) = T(end_ind)-T(start_ind);
% end
% 
% t_ampa = 0.005;
% t_nmda = 0.1;
% tau_ie = 0.5*t_ampa+0.5*t_nmda;
% tau_ee = (0.5-qdiff)*t_ampa+(0.5+qdiff)*t_nmda;
% delta_tau = (tau_ee - tau_ie)*10^3;
% 
% figure('color','w')
% plot(stim_hold,decay_all,'--k','LineWidth',2)
% hold on
% plot(stim_hold,decay_qd_all,'k','LineWidth',2)
% set(gca,'FontSize',20)
% legend('p=0.1',['\Delta\tau=',num2str(delta_tau),' ms'])
% xlabel('Stimulus length (s)','FontSize',30)
% ylabel('Decay time (s)','FontSize',30)
% 
% 
% %%%% Plot decay time as a function of steady state firing rate
% stim = 2.0;
% isi = 2.0;
% num_reps = 1;
% stim_setup = [stim,isi,num_reps];
% 
% tr_ff = 0.5;
% U_ff = 0.2;
% 
% qdiff = 0;
% 
% percent = 0.1;
% delta_U = percent*U_base;
% delta_tr = percent*tr_base;
% 
% U_base = 0.1;
% U = U_base*ones(size(q));
% U = U+delta_U*diff_template;
% 
% tr_base = 0.5;
% tr = tr_base*ones(size(q));
% tr = tr+delta_tr*diff_template;
% 
% all_hz = 1:1:50;
% p_win = zeros(length(all_hz),1);
% for i = 1:length(all_hz)
%     r = all_hz(i);
%     win_ss = @(w_in)(dfdbk_ffstd_steadystate(w,k,tr,U,tr_ff,U_ff,I_step,w_in)-r);
%     p_win(i) = fzero(win_ss,[0.01,100]);
% end
% 
% decay_p_all = zeros(length(all_hz),1);
% 
% parfor i = 1:length(all_hz)
%     w_in = p_win(i);
%     [rise_time, decay_time] = dfdbk_ffstd_manysynapses_RandD(w,k,q,qdiff,tr,U,tr_ff,U_ff,I_step,w_in,q_in,stim_setup);
%     decay_p_all(i) = decay_time;
% end
% 
% 
% 
% 
% solutions = cell(3,1);
% for i = 1:3
%     w_in = p_win(i*10);
%     [sol, ss_all] = dfdbk_ffstd_manysynapses_sim(w,k,q,qdiff,tr,U,tr_ff,U_ff,I_step,w_in,q_in,stim_setup);
%     solutions{i} = sol;
% end
% 
% figure('color','w')
% for i = 1:3
%     holder = solutions{i};
%     T = holder.x;
%     ye = holder.y(1,:);
%     plot(T,ye,'LineWidth',2)
%     hold on
% end
% set(gca,'FontSize',20)
% legend(['R_e^{ss}=10 Hz'],['R_e^{ss}=20 Hz'],['R_e^{ss}=30 Hz'])
% xlabel('Time (s)','FontSize',30)
% ylabel('Activity (Hz)','FontSize',30)
% 
% 
% qdiff = 0.05;
% 
% percent = 0.0;
% delta_U = percent*U_base;
% delta_tr = percent*tr_base;
% 
% U_base = 0.1;
% U = U_base*ones(size(q));
% U = U+delta_U*diff_template;
% 
% tr_base = 0.5;
% tr = tr_base*ones(size(q));
% tr = tr+delta_tr*diff_template;
% 
% qdiff_win = zeros(length(all_hz),1);
% for i = 1:length(all_hz)
%     r = all_hz(i);
%     win_ss = @(w_in)(dfdbk_ffstd_steadystate(w,k,tr,U,tr_ff,U_ff,I_step,w_in)-r);
%     qdiff_win(i) = fzero(win_ss,[0.01,100]);
% end
% 
% decay_qdiff_all = zeros(length(all_hz),1);
% 
% parfor i = 1:length(all_hz)
%     w_in = qdiff_win(i);
%     [rise_time, decay_time] = dfdbk_ffstd_manysynapses_RandD(w,k,q,qdiff,tr,U,tr_ff,U_ff,I_step,w_in,q_in,stim_setup);
%     decay_qdiff_all(i) = decay_time;
% end
% 
% t_ampa = 0.005;
% t_nmda = 0.1;
% tau_ie = 0.5*t_ampa+0.5*t_nmda;
% tau_ee = (0.5-qdiff)*t_ampa+(0.5+qdiff)*t_nmda;
% delta_tau = (tau_ee - tau_ie)*10^3;
% 
% figure('color','w')
% plot(all_hz,decay_p_all,'--k','LineWidth',2)
% hold on
% plot(all_hz,decay_qdiff_all,'k','LineWidth',2)
% set(gca,'FontSize',20)
% legend('p=0.1',['\Delta\tau=',num2str(delta_tau),' ms'])
% xlabel('R_e^{ss}','FontSize',30)
% ylabel('Decay time (s)','FontSize',30)

% solutions = cell(3,1);
% for i = 1:3
%     w_in = qdiff_win(i*10);
%     [sol, ss_all] = dfdbk_ffstd_manysynapses_sim(w,k,q,qdiff,tr,U,tr_ff,U_ff,I_step,w_in,q_in,stim_setup);
%     solutions{i} = sol;
% end
% 
% figure('color','w')
% for i = 1:3
%     holder = solutions{i};
%     T = holder.x;
%     ye = holder.y(1,:);
%     plot(T,ye,'LineWidth',2)
%     hold on
% end
% set(gca,'FontSize',14)
% legend(stim_leg{:})
% xlabel('Time (s)','FontSize',25)
% ylabel('Activity (Hz)','FontSize',25)
% title('qdiff')





% [sol, ss_all] = dfdbk_std_manysynapses_sim(w,k,q,qdiff,tr,U,I_step,w_in,q_in,stim_setup);
% [rise_time, decay_time] = dfdbk_std_manysynapses_RandD(w,k,q,qdiff,tr,U,I_step,w_in,q_in,stim_setup)
% 
% 
% % cmp = get(gca,'colororder');
% 
% 
% for i = 1:1
%     figure('color','w')
%     plot(sol.x,sol.y(i,:),'LineWidth',2)
%     hold on
%     plot([sol.x(1),sol.x(end)],[ss_all(i),ss_all(i)],'--k')
% end
% set(gca,'FontSize',14)
% xlabel('Time (s)','FontSize',25)
% ylabel('Activity (Hz)','FontSize',25)


%%%%%%%%%%%%%%%%
% Plot delta tau and compute associated u_base
%%%%%%%%%%%%%%%%


% u = 0.01:0.001:0.4;
% percent = 0.05;
% tr_ss = 0.5;
% R = 20;
% x = @(u,percent)1./(1+(1+percent)^2*u*tr_ss*R)
% q = [0.25,0.75]'
% f1 = @(u_base,percent)(q(2)-q(1))*(0.005-0.1)*(x(u_base,percent)-x(u_base,-percent))./(x(u_base,percent)+x(u_base,-percent));
% 
% figure('color','w')
% plot(u,1000*f1(u,percent))
% ylabel('\Delta\tau (ms)')
% xlabel('u')
% 
% f2 = @(u_base,percent)(q(2)-q(1))*(0.005-0.1)*(x(u_base,percent)-x(u_base,-percent));
% 
% figure('color','w')
% plot(u,100*f2(u,percent))
% ylabel('\Delta\tau (ms)')
% xlabel('u')
% 
% combined_func = f2(u,percent);
% [max_val, mid_ind] = max(combined_func);
% mid_u = u(mid_ind);
% [na, small_ind] = find(combined_func>((9/10)*max_val),1,'first');
% small_u = u(small_ind);
% [na, large_ind] = find(combined_func>((9/10)*max_val),1,'last');
% large_u = u(large_ind);
% U_base = {small_u,mid_u,large_u};
% 
% 
% tr_ff = 0.5;
% U_ff = 0.15;
% 
% U_leg = {['u=',num2str(U_base{1})],['u=',num2str(U_base{2})],['u=',num2str(U_base{3})]}
% solutions = cell(length(U_base),1);
% for i = 1:length(U_base)
%     U = U_base{i}*ones(size(q));
%     U = U+delta_U*diff_template;
%     win_ss = @(w_in)(dfdbk_ffstd_steadystate(w,k,tr,U,tr_ff,U_ff,I_step,w_in)-20.0);
%     w_in = fzero(win_ss,[0.1,10]);
% 
%     [sol, ss_all] = dfdbk_ffstd_manysynapses_sim(w,k,q,qdiff,tr,U,tr_ff,U_ff,I_step,w_in,q_in,stim_setup);
%     [rise_time, decay_time] = dfdbk_ffstd_manysynapses_RandD(w,k,q,qdiff,tr,U,tr_ff,U_ff,I_step,w_in,q_in,stim_setup)
%     solutions{i} = sol;
% end
% 
% % cmp = get(gca,'colororder');
% 
% figure('color','w')
% for i = 1:length(U_base)
%     plot(solutions{i}.x,solutions{i}.y(1,:),'LineWidth',2)
%     hold on
% end
% set(gca,'FontSize',14)
% xlabel('Time (s)','FontSize',25)
% ylabel('Activity (Hz)','FontSize',25)
% legend(U_leg{:})

% for i = 1:1
%     figure('color','w')
%     plot(sol.x,sol.y(i,:),'LineWidth',2)
%     hold on
%     plot([sol.x(1),sol.x(end)],[ss_all(i),ss_all(i)],'--k')
% end
% set(gca,'FontSize',14)
% xlabel('Time (s)','FontSize',25)
% ylabel('Activity (Hz)','FontSize',25)



























