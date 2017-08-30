% Produces the sensitivity analysis plots.

% I_step = 15.0;
% 
% w = 100;
% k = 1.1;
% 
% w_in = 1.0;
% q_in = 0.5;
% 
% 
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
% msize=16;
% 
% %%%% Rise and decay across many percent and delta_q values.
% 
% percent_hold = [0:0.005:0.15];
% delta_q = [0.0,0.1,0.25,0.4];
% rise_all = zeros(length(percent_hold),length(delta_q));
% decay_all = zeros(length(percent_hold),length(delta_q));
% 
% for j = 1:length(delta_q)
%     q = [0.5-delta_q(j),0.5+delta_q(j)]';
%     parfor i = 1:length(percent_hold)
%         percent = percent_hold(i);
%         delta_U = percent*U_base;
%         delta_tr = percent*tr_base;
% 
%         U = U_base*ones(size(q));
%         U = U+delta_U*diff_template;
% 
%         tr = tr_base*ones(size(q));
%         tr = tr+delta_tr*diff_template;
% 
%         [rise_time, decay_time] = dfdbk_std_manysynapses_RandD(w,k,q,qdiff,tr,U,I_step,w_in,q_in,stim_setup);
%         rise_all(i,j) = rise_time;
%         decay_all(i,j) = decay_time;
%     end
% end
% 
% figure('color','w')
% for i = 1: length(delta_q)
%     plot(percent_hold,rise_all(:,i),'LineWidth',2)
%     hold on
% end
% set(gca,'FontSize',20)
% legend('\Delta q= 0.0','\Delta q= 0.1','\Delta q= 0.25','\Delta q= 0.4','location','northwest')
% 
% ind1 = find(percent_hold>=0.07,1);
% ind2 = find(percent_hold>=0.05,1);
% hLine = plot(0.07,rise_all(ind1,3),'k+','markersize',msize,'LineWidth',2);
% set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% hold on
% hLine = plot(0.05,rise_all(ind2,3),'kx','markersize',msize,'LineWidth',2);
% set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% 
% ylim([0,0.06])
% xlabel('p','FontSize',30)
% ylabel('Rise Time (s)','FontSize',30)
% 
% 
% figure('color','w')
% for i = 1:length(delta_q)
%     plot(percent_hold,decay_all(:,i),'LineWidth',2)
%     hold on
% end
% set(gca,'FontSize',20)
% legend('\Delta q= 0.0','\Delta q= 0.1','\Delta q= 0.25','\Delta q= 0.4','location','northwest')
% 
% ind1 = find(percent_hold>=0.07,1);
% ind2 = find(percent_hold>=0.05,1);
% hLine = plot(0.07,decay_all(ind1,3),'k+','markersize',msize,'LineWidth',2);
% set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% hold on
% hLine = plot(0.05,decay_all(ind2,3),'kx','markersize',msize,'LineWidth',2);
% set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% 
% xlabel('p','FontSize',30)
% ylabel('Decay Time (s)','FontSize',30)
% 
% 
% 
% %%%% Rise and decay across many percent and q_shift values.
% 
% percent_hold = [0:0.005:0.15];
% q_shift = [-0.2,-0.1,0.0,0.1,0.2];
% rise_all = zeros(length(percent_hold),length(q_shift));
% decay_all = zeros(length(percent_hold),length(q_shift));
% 
% for j = 1:length(q_shift)
%     q = [0.25+q_shift(j),0.75+q_shift(j)]';
%     parfor i = 1:length(percent_hold)
%         percent = percent_hold(i);
%         delta_U = percent*U_base;
%         delta_tr = percent*tr_base;
% 
%         U = U_base*ones(size(q));
%         U = U+delta_U*diff_template;
% 
%         tr = tr_base*ones(size(q));
%         tr = tr+delta_tr*diff_template;
% 
%         [rise_time, decay_time] = dfdbk_std_manysynapses_RandD(w,k,q,qdiff,tr,U,I_step,w_in,q_in,stim_setup);
%         rise_all(i,j) = rise_time;
%         decay_all(i,j) = decay_time;
%     end
% end
% 
% figure('color','w')
% for i = 1:length(q_shift)
%     plot(percent_hold,rise_all(:,i),'LineWidth',2)
%     hold on
% end
% set(gca,'FontSize',20)
% legend('q_{shift}= -0.2','q_{shift}= -0.1','q_{shift}= 0.0','q_{shift}= 0.1','q_{shift}= 0.2','location','northwest')
% 
% 
% hLine = plot(0.07,rise_all(ind1,3),'k+','markersize',msize,'LineWidth',2);
% set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% hold on
% hLine = plot(0.05,rise_all(ind2,3),'kx','markersize',msize,'LineWidth',2);
% set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% 
% ylim([0,0.06])
% xlabel('p','FontSize',30)
% ylabel('Rise Time (s)','FontSize',30)
% 
% figure('color','w')
% for i = 1:length(q_shift)
%     plot(percent_hold,decay_all(:,i),'LineWidth',2)
%     hold on
% end
% set(gca,'FontSize',20)
% legend('q_{shift}= -0.2','q_{shift}= -0.1','q_{shift}= 0.0','q_{shift}= 0.1','q_{shift}= 0.2','location','northwest')
% 
% ind1 = find(percent_hold>=0.07,1);
% ind2 = find(percent_hold>=0.05,1);
% hLine = plot(0.07,decay_all(ind1,3),'k+','markersize',msize,'LineWidth',2);
% set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% hold on
% hLine = plot(0.05,decay_all(ind2,3),'kx','markersize',msize,'LineWidth',2);
% set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% 
% xlabel('p','FontSize',30)
% ylabel('Decay Time (s)','FontSize',30)




%%%% Changing p values.

I_step = 15.0;

w = 100;
k = 1.1;

w_in = 1.0;
q_in = 0.5;


U_base = 0.1;
tr_base = 0.5;
qdiff = -0.0075;
diff_template = [ones(1,1);-ones(1,1)];

stim = 2.0;
isi = 2.0;
num_reps = 1;
stim_setup = [stim,isi,num_reps];

msize=16;

%%%% Rise and decay across many percent and delta_q values.

percent1_hold = [0:0.005:0.15];
percent2_hold = [0:0.005:0.15];
rise_all = zeros(length(percent_hold),length(delta_q));
decay_all = zeros(length(percent_hold),length(delta_q));

parfor i = 1:length(percent_hold)
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

    [rise_time, decay_time] = dfdbk_std_manysynapses_changep_RandD(w,k,q,qdiff,tr,U,I_step,w_in,q_in,stim_setup);
    rise_all(i) = rise_time;
    decay_all(i) = decay_time;
end

figure('color','w')
for i = 1: length(delta_q)
    plot(percent_hold,rise_all(:,i),'LineWidth',2)
    hold on
end
set(gca,'FontSize',20)
legend('\Delta q= 0.0','\Delta q= 0.1','\Delta q= 0.25','\Delta q= 0.4','location','northwest')

ind1 = find(percent_hold>=0.07,1);
ind2 = find(percent_hold>=0.05,1);
hLine = plot(0.07,rise_all(ind1,3),'k+','markersize',msize,'LineWidth',2);
set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold on
hLine = plot(0.05,rise_all(ind2,3),'kx','markersize',msize,'LineWidth',2);
set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

ylim([0,0.06])
xlabel('p','FontSize',30)
ylabel('Rise Time (s)','FontSize',30)


