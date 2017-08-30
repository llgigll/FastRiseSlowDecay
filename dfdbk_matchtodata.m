% Setup for monkey 1

I_step = 15.0;

w = 100;
k = 1.1;

w_in = 1.0;
q_in = 0.5;


q = [0.25,0.75]';
U_base = 0.1;
tr_base = 0.5;
diff_template = [ones(1,1);-ones(1,1)];

latency_offset = 0.103-0.055;
stim = 0.5+latency_offset;
isi = 1.0-latency_offset;
num_reps = 1;
stim_setup = [stim,isi,num_reps];

tr_ff = 0.5;
U_ff = 0.1;

qdiff = -0.004;

%percent = 0.07;
percent = 0.05;
delta_U = percent*U_base;
delta_tr = percent*tr_base;

U_base = 0.1;
U = U_base*ones(size(q));
U = U+delta_U*diff_template;

tr_base = 0.5;
tr = tr_base*ones(size(q));
tr = tr+delta_tr*diff_template;

all_hz = 1:1:20;
p_win = zeros(length(all_hz),1);
for i = 1:length(all_hz)
    r = all_hz(i);
    win_ss = @(w_in)(dfdbk_ffstd_steadystate(w,k,tr,U,tr_ff,U_ff,I_step,w_in)-r);
    p_win(i) = fzero(win_ss,[0.01,100]);
end

decay_p_all = zeros(length(all_hz),1);

parfor i = 1:length(all_hz)
    w_in = p_win(i);
    [rise_time, decay_time] = dfdbk_ffstd_manysynapses_RandD(w,k,q,qdiff,tr,U,tr_ff,U_ff,I_step,w_in,q_in,stim_setup);
    decay_p_all(i) = decay_time;
end




solutions = cell(1,1);
for i = 1:1
    w_in = p_win(8);
    [sol, ss_all] = dfdbk_ffstd_manysynapses_sim(w,k,q,qdiff,tr,U,tr_ff,U_ff,I_step,w_in,q_in,stim_setup);
    solutions{i} = sol;
end

load AE2papercells.mat
holder = origpop.AmbigEdge2.analysis.difference.allcellsrates(:,:,1:65);
holder = mean(holder,3);
for i = 9:9
    Pe_rate = holder(:,i);
    filter = fspecial('gaussian',[30 1],6.0); % gaussian kernel where s= size of contour
    Pe_rate = conv(Pe_rate, filter,'same');
    figure('color','w')
    plot((301:2:2001)/1000-0.555,Pe_rate(151:1001))
    hold on
    xlim(([555,2001]/1000-0.555))
end
hold on

for i = 1:1
    holder = solutions{i};
    T = holder.x;
    ye = holder.y(1,:);
    plot(T,ye,'LineWidth',2)
    hold on
end
set(gca,'FontSize',20)
legend(['Monkey TH'],['Model'])
xlabel('Time (s)','FontSize',30)
ylabel('Activity (Hz)','FontSize',30)