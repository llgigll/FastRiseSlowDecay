% Match data to different input lengths

load AEFlashAE2combo.mat
holder = test1pop.AmbigEdgeFlash.analysis.difference.ratesbinbybin;
% holder = mean(holder,3);
figure('color','w')
ind = [1,4];
start = [-0.082,0.038];
for j = 1:2
    i = ind(j); 
    Pe_rate = holder(:,i);
    filter = fspecial('gaussian',[30 1],8.0); % gaussian kernel where s= size of contour
    Pe_rate = conv(Pe_rate, filter,'same');
    plot((1:2:3001)/1000+start(j)-0.545,Pe_rate(1:1501))
    hold on
end


holder = test2pop.AmbigEdge2.analysis.difference.ratesbinbybin;
% holder = mean(holder,3);

ind = [1:9];
% start = [-150,-50]
for j = 9:9
    i = ind(j); 
    Pe_rate = holder(:,i);
    filter = fspecial('gaussian',[30 1],8.0); % gaussian kernel where s= size of contour
    Pe_rate = conv(Pe_rate, filter,'same');
    plot((1:2:3001)/1000-0.58,Pe_rate(1:1501))
    hold on
    plot([1,3001]/1000-1.045,[0,0],'--k')
    xlim(([1,3001]/1000-0.5))
end


I_step = 15.0;

w = 100;
k = 1.1;
q_in = 0.5;

q = [0.25,0.75]';
tr_base = 0.5;
diff_template = [ones(1,1);-ones(1,1)];

tr_ff = 0.5;
U_ff = 0.1;
U_base = 0.1;
qdiff = -0.004;
w_in_hold = [0.3,0.25,0.3]*(1+tr_ff*U_ff*I_step);

percent = 0.1;
delta_U = percent*U_base;
delta_tr = percent*tr_base;

U = U_base*ones(size(q));
U = U+delta_U*diff_template;

tr = tr_base*ones(size(q));
tr = tr+delta_tr*diff_template;


stim_hold = [0.25,0.15,0.5];
isi = 2.0;
num_reps = 1;

solutions = cell(length(stim_hold),1);
for i = 1:length(stim_hold)
    w_in = w_in_hold(i);
    stim_setup = [stim_hold(i),isi,num_reps];
    [sol, ss_all] = dfdbk_ffstd_manysynapses_sim(w,k,q,qdiff,tr,U,tr_ff,U_ff,I_step,w_in,q_in,stim_setup);
    solutions{i} = sol;
end

percent_leg = cell(length(stim_hold),1);
for i = 1:length(stim_hold)
   stim_leg{i} = ['Stimulus Length=',num2str(stim_hold(i)),' s'];
end

figure('color','w')
for i = 1:length(stim_hold)
    holder = solutions{i};
    T = holder.x;
    ye = holder.y(1,:);
    plot(T-stim_hold(i)+0.5,ye,'LineWidth',2)
    hold on
end
set(gca,'FontSize',20)
legend(stim_leg{:})
xlabel('Time (s)','FontSize',30)
ylabel('Activity (Hz)','FontSize',30)