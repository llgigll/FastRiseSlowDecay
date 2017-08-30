%% Load spiking files from Brian
%  Loads a set of .mat files produced in Brian and averages and plots them.
    
figure('Color','w');
u = {'20','10','5','0'};
u_leg = {'u=0.20','u=0.10','u=0.05','No STD'};
for i = 1:length(u) 

    load(['./python/pfdbk_final4_u',u{i},'_network3.mat'])
    
    filter = fspecial('gaussian',[30 1],4.0); % gaussian kernel where s= size of contour

    Pe_rate = conv(Pe_rate, filter,'same');
    plot(Pe_time-0.5,Pe_rate,'LineWidth',2)
    hold on

end
xlim([0,4])
set(gca,'fontsize',20)
xlabel('Time (s)','FontSize',25)
ylabel('Rate (Hz)','FontSize',25)
legend(u_leg{:})
% 
% 
% 
% %  Plot impact of STD on dfdbk
% 
% figure('Color','w');
% u = {'20','10','5','0'};
% u_leg = {'u=0.20','u=0.10','u=0.05','No STD'};
% for i = 1:length(u) 
% 
%     load(['./python/dfdbk_final_u',u{i},'_network3.mat'])
%     
%     filter = fspecial('gaussian',[30 1],4.0); % gaussian kernel where s= size of contour
% 
%     Pe_rate = conv(Pe_rate, filter,'same');
%     plot(Pe_time-0.5,Pe_rate,'LineWidth',2)
%     hold on
% 
% end
% xlim([0,4])
% xlabel('Time (s)','FontSize',30)
% ylabel('Rate (Hz)','FontSize',30)
% set(gca,'fontsize',20)
% legend(u_leg{:})

%  Plot increasing delta q

figure('Color','w');
delta_q = {'0','10','20'};
diff_plot = [0,0.1,0.2];
t_ampa = 0.01;
t_nmda = 0.2;
dtau = zeros(length(diff_plot),1);
for i = 1:length(dtau)
    tau_ie = (0.5+diff_plot(i))*t_ampa+(0.5-diff_plot(i))*t_nmda;
    tau_ee = (0.5-diff_plot(i))*t_ampa+(0.5+diff_plot(i))*t_nmda;
    dtau(i) = (tau_ee - tau_ie)*10^3;
end
delta_tau_leg = {['\Delta\tau=',num2str(dtau(1)),' ms'],['\Delta\tau=',num2str(dtau(2)),' ms'],['\Delta\tau=',num2str(dtau(3)),' ms']};
for i = 1:length(delta_q) 

    load(['./python/dfdbk_final2_deltaq',delta_q{i},'_network3.mat'])
    
    filter = fspecial('gaussian',[30 1],4.0); % gaussian kernel where s= size of contour

    Pe_rate = conv(Pe_rate, filter,'same');
    plot(Pe_time-0.5,Pe_rate,'LineWidth',2)
    hold on

end
xlim([-0.1,4])
ylim([0,45])
xlabel('Time (s)','FontSize',30)
ylabel('Rate (Hz)','FontSize',30)
set(gca,'fontsize',20)
legend(delta_tau_leg{:})




%  Plot increasing delta q with ffd std

figure('Color','w');
uff = {'10','25','50'};
uff_leg = {'u_{ff}=0.10','u_{ff}=0.25','u_{ff}=0.50'};
for i = 1:length(uff) 

    load(['./python/dfdbk_final2_ffdstd_uff_',uff{i},'_network3.mat'])
    
    filter = fspecial('gaussian',[30 1],4.0); % gaussian kernel where s= size of contour

    Pe_rate = conv(Pe_rate, filter,'same');
    plot(Pe_time-0.5,Pe_rate,'LineWidth',2)
    hold on

end
xlim([-0.1,4])
ylim([0,45])
xlabel('Time (s)','FontSize',30)
ylabel('Rate (Hz)','FontSize',30)
set(gca,'fontsize',20)
legend(uff_leg{:})




%  Plot increasing p

figure('Color','w');
p = {'0','12','25'};
p_leg = {'p=0','p=0.125','p=0.25'};
for i = 1:length(p) 

    load(['./python/dfdbk_final2_p_',p{i},'_network3.mat'])
    
    filter = fspecial('gaussian',[30 1],4.0); % gaussian kernel where s= size of contour

    Pe_rate = conv(Pe_rate, filter,'same');
    plot(Pe_time-0.5,Pe_rate,'LineWidth',2)
    hold on

end
xlim([-0.1,4])
ylim([0,45])
xlabel('Time (s)','FontSize',30)
ylabel('Rate (Hz)','FontSize',30)
set(gca,'fontsize',20)
legend(p_leg{:})






%  Plot increasing p with ffd std

% figure('Color','w');
% p = {'0','12','25'};
% for i = 1:length(p) 
% 
%     load(['./python/dfdbk_final_ffdstd_p_',p{i},'_network3.mat'])
%     
%     filter = fspecial('gaussian',[30 1],4.0); % gaussian kernel where s= size of contour
% 
%     Pe_rate = conv(Pe_rate, filter,'same');
%     plot(Pe_time-0.5,Pe_rate,'LineWidth',2)
%     hold on
% 
% end
% xlim([0,4])
% xlabel('Time (s)','FontSize',30)
% ylabel('Rate (Hz)','FontSize',30)
% set(gca,'fontsize',20)
% legend(p_leg{:})


%  Plot increasing p with ffd std

figure('Color','w');
uff = {'10','25','50'};
uff_leg = {'u_{ff}=0.10','u_{ff}=0.25','u_{ff}=0.50'};
for i = 1:length(uff) 

    load(['./python/dfdbk_final_ffdstd_p0125_uff',uff{i},'_network3.mat'])
    
    filter = fspecial('gaussian',[30 1],4.0); % gaussian kernel where s= size of contour

    Pe_rate = conv(Pe_rate, filter,'same');
    plot(Pe_time-0.5,Pe_rate,'LineWidth',2)
    hold on

end
xlim([-0.1,4])
ylim([0,45])
xlabel('Time (s)','FontSize',30)
ylabel('Rate (Hz)','FontSize',30)
set(gca,'fontsize',20)
legend(uff_leg{:})


% Testing time constant Plot impact of STD on dfdbk
% 
% figure('Color','w');
% u = {'5'};
% u_leg = {'orig figure','tau200','tau100',''};
% 
% 
% load(['./python/dfdbk_final_u5_network3.mat'])
% 
% filter = fspecial('gaussian',[30 1],4.0); % gaussian kernel where s= size of contour
% 
% Pe_rate = conv(Pe_rate, filter,'same');
% plot(Pe_time-0.5,Pe_rate,'LineWidth',2)
% hold on
% 
% load(['./python/dfdbk_testing_u5_network3.mat'])
% 
% filter = fspecial('gaussian',[30 1],4.0); % gaussian kernel where s= size of contour
% 
% Pe_rate = conv(Pe_rate, filter,'same');
% plot(Pe_time-0.5,Pe_rate,'LineWidth',2)
% hold on
% 
% load(['./python/dfdbk_testing2_u5_network3.mat'])
% 
% filter = fspecial('gaussian',[30 1],4.0); % gaussian kernel where s= size of contour
% 
% Pe_rate = conv(Pe_rate, filter,'same');
% plot(Pe_time-0.5,Pe_rate,'LineWidth',2)
% 
% 
% xlim([0,4])
% xlabel('Time (s)','FontSize',30)
% ylabel('Rate (Hz)','FontSize',30)
% set(gca,'fontsize',20)
% legend(u_leg{:})





diff_plot = [-0.0065,-0.0085,-0.011,-0.003];
t_ampa = 0.01;
t_nmda = 0.2;
dtau = zeros(length(diff_plot),1);
for i = 1:length(dtau)
    tau_ie = (0.5+diff_plot(i))*t_ampa+(0.5-diff_plot(i))*t_nmda;
    tau_ee = (0.5-diff_plot(i))*t_ampa+(0.5+diff_plot(i))*t_nmda;
    dtau(i) = (tau_ee - tau_ie)*10^3;
end















