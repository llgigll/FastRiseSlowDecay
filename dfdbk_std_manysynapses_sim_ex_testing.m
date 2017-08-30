% Produces plots for figures 6 and 8.

% Plot three trajectories of the derivative feedback network with STD on
% the excitatory to excitatory projections.

% I_step = 15.0;
% 
% w = 100;
% k = 1.1;
% 
% w_in = 1.0;
% q_in = 0.5;
% 
% 
% q = [0.25,0.75]';
% U_base = 0.1;
% tr_base = 0.5;
% qdiff = -0.0075;
% diff_template = [ones(1,1);-ones(1,1)];
% 
% stim = 2.0;
% isi = 2.0;
% num_reps = 1;
% stim_setup = [stim,isi,num_reps];
% 
% %%%%  Simulations for three different percent values.
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
%     
%     [sol, ss_all] = dfdbk_std_manysynapses_sim(w,k,q,qdiff,tr,U,I_step,w_in,q_in,stim_setup);
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
% figure('color','w')
% for i = 1:length(percent_hold)
%     holder = solutions{i};
%     x1 = holder.y(end-1,:);
%     x2 = holder.y(end,:);
%     delta_tau = (q(2)-q(1))*(0.005-0.1)*(x1-x2)./(x1+x2)+qdiff*(0.1-0.005);
%     T = holder.x;
%     plot(T,delta_tau*1000,'LineWidth',2)
%     hold on
% end
% set(gca,'FontSize',20)
% legend(percent_leg{:})
% xlabel('Time (s)','FontSize',30)
% ylabel('\Delta\tau (ms)','FontSize',30)


%%%%%%%%%%%%%%%%%%%%% This allows for 2 different values of p %%%%%%%%%%%%%

I_step = 15.0;

w = 100;
k = 1.1;

w_in = 1.0;
q_in = 0.5;


q = [0.25,0.75]';
U_base = 0.1;
tr_base = 0.5;
% qdiff = -0.0075;
qdiff = -0.0075;
diff_template = [ones(1,1);-ones(1,1)];

stim = 2.0;
isi = 2.0;
num_reps = 1;
stim_setup = [stim,isi,num_reps];

%%%%  Simulations for three different percent values.
percent1_hold = [0.0,0.05,0.1,0.15]*1.0;
percent2_hold = -[0.0,0.05,0.1,0.15]*0.0;

% percent2_hold = -[0.0,0.00,0.0,0.0];
solutions = cell(length(percent_hold),1);
for i = 1:length(percent1_hold)
    percent1 = percent1_hold(i);
    delta_U = percent1*U_base;
    delta_tr = percent1*tr_base;
    U1 = U_base*ones(size(q));
    U1 = U1+delta_U*diff_template;
    tr1 = tr_base*ones(size(q));
    tr1 = tr1+delta_tr*diff_template;
    
    percent2 = percent2_hold(i);
    delta_U = percent2*U_base;
    delta_tr = percent2*tr_base;
    U2 = U_base*ones(size(q));
    U2 = U2+delta_U*diff_template;
    tr2 = tr_base*ones(size(q));
    tr2 = tr2+delta_tr*diff_template;
    
    U = [U1;U2];
    tr = [tr1;tr2];
    
    [sol, ss_all] = dfdbk_std_manysynapses_changep_sim(w,k,q,qdiff,tr,U,I_step,w_in,q_in,stim_setup);
    solutions{i} = sol;
end

percent_leg = cell(length(percent1_hold),1);
for i = 1:length(percent1_hold)
   percent_leg{i} = ['p=',num2str(percent1_hold(i))];
end

figure('color','w')
for i = 1:length(percent1_hold)
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




%%%%%%%%%%%%%%%%%%%%% This allows for 2 different values of p %%%%%%%%%%%%%

I_step = 15.0;

w = 100;
k = 1.1;

w_in = 1.0;
q_in = 0.5;


q = [0.25,0.75]';
U_base = 0.1;
tr_base = 0.5;
% qdiff = -0.0075;
qdiff = -0.0075;
diff_template = [ones(1,1);-ones(1,1)];

stim = 2.0;
isi = 2.0;
num_reps = 1;
stim_setup = [stim,isi,num_reps];

%%%%  Simulations for three different percent values.
percent1_hold = [0.0,0.05,0.1,0.15]*1.0;
R_ss=25;
% percent2_hold = -[0.0,0.00,0.0,0.0];
solutions = cell(length(percent_hold),1);
for i = 1:length(percent1_hold)
    percent1 = percent1_hold(i);
    delta_U = percent1*U_base;
    delta_tr = percent1*tr_base;
    U1 = U_base*ones(size(q));
    U1 = U1+delta_U*diff_template;
    tr1 = tr_base*ones(size(q));
    tr1 = tr1+delta_tr*diff_template;
    
    x1 = 1/(1+U1(1)*tr1(1)*R_ss);
    x2 = 1/(1+U1(2)*tr1(2)*R_ss);
    x_ave = (x1+x2)/2;
    
    pie = sqrt((1-x_ave)/(x_ave*U_base*tr_base*R_ss))-1
    
    U2 = [U_base*(1+pie);U_base*(1+pie)]*1.00;
    tr2 = [tr_base*(1+pie);tr_base*(1+pie)]*1.00;
    
    U = [U1;U2];
    tr = [tr1;tr2];
    
    [sol, ss_all] = dfdbk_std_manysynapses_changep_sim(w,k,q,qdiff,tr,U,I_step,w_in,q_in,stim_setup);
    solutions{i} = sol;
end

percent_leg = cell(length(percent1_hold),1);
for i = 1:length(percent1_hold)
   percent_leg{i} = ['p=',num2str(percent1_hold(i))];
end

figure('color','w')
for i = 1:length(percent1_hold)
    holder = solutions{i};
    T = holder.x;
    ye = holder.y(1,:);
    T_ind = find(T>=0.99*stim,1);
    ye_hold = ye(T_ind);
    ye_hold = 1;
    plot(T,ye/ye_hold,'LineWidth',2)
    hold on
end
set(gca,'FontSize',20)
legend(percent_leg{:})
xlabel('Time (s)','FontSize',30)
ylabel('Activity (Hz)','FontSize',30)

