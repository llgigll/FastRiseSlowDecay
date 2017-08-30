function sol = dfdbk_amplification_risetime()



Tfinal = 2.5;
Ie = 1;
Ii = 0.0;
k = [1.05,1.1,1.25,1.5,2.0];
k_leg =  cell(length(k),1);
for i = 1:length(k)
    k_leg{i} = ['k = ',num2str(k(i))];
end
w = 1:1:60;
% f = @(k,w)(1+k*w)./(1+(k-1)*w);
f = @(k,w,Ie,Ii) (Ie*(1+k*w)-Ii*k*w)./(1+(k-1)*w);

figure('color','w')
for i = k; plot(w,f(i,w,Ie,Ii),'LineWidth',2), hold on, end
ylabel('Gain','FontSize',30)
xlabel('w','FontSize',30)
ax = gca; 
ax.FontSize = 16;
legend(k_leg{:},'location','northwest')
clear i

diff = [-0.01,0.0,0.05];
diff_leg =  cell(length(diff),1);
for i = 1:length(diff)
    diff_leg{i} = ['\Delta','q = ',num2str(diff(i))];
end
%diff = 0.0;
%diff = -0.01;
q = 0.5;
for m = 1:length(diff)
    rise = zeros(length(k),length(w));
    for i = 1:length(k)
        for j = 1:length(w)
            rise(i,j) = risetime(diff(m), q, w(j), k(i));
        end
    end
    
    figure('color','w')
    for i = 1:length(k)
        plot(w,rise(i,:),'LineWidth',2)
        hold on
    end
    ylabel('Rise time (s)','FontSize',30)
    xlabel('w','FontSize',30)
    ax = gca; 
    ax.FontSize = 16;
end
clear i j m

k_plot = 1.25;
w_plot = 30;
figure('color','w')
lines = {'k-','k-.','k--'};
for i = 1:length(diff)
    sys = system(diff(i),q,w_plot,k_plot);
    [y,t] = step(sys,Tfinal);
    plot(t,y,lines{i},'LineWidth',2)
    hold on
end
ylabel('Time (s)','FontSize',30)
xlabel('Activity (Hz)','FontSize',30)
ax = gca; 
ax.FontSize = 16;
legend(diff_leg{:},'location','southeast')

k_plot = 1.25;
w_plot = 30;
for i = 1:length(diff)
    figure('color','w')
    sys = system(diff(i),q,w_plot,k_plot);
    bode(sys)
end

    
% 
% figure('color','w')
% for i = 1:length(k)
%     plot(w,f(k(i),w)./rise(i,:))
%     hold on
% end
% ylabel('Gain / Rise time (s)')
% xlabel('w')



 function rise = risetime(diff,q,w,k)
    te  = 0.02;
    ti  = 0.01;
    tee_a = 0.005;
    tee_n = 0.1;
    tie_a = 0.005;
    tie_n = 0.1;
    tei = 0.01;
    tii = 0.01;
    qdiff1 = diff;
    qdiff2 = 0;
    Jee = w;
    Jei = k*w;
    Jie = w;
    Jii = k*w;
    fpos = 0;

    A = [-1/te,0,((1-q)-qdiff1)*(Jee+fpos)/te,(q+qdiff1)*(Jee+fpos)/te,0,0,-Jei/te,0;...
         0,-1/ti,0,0,((1-q)-qdiff2)*Jie/ti,(q+qdiff2)*Jie/ti,0,-Jii/ti;...
         1/tee_a,0,-1/tee_a,0,0,0,0,0;...
         1/tee_n,0,0,-1/tee_n,0,0,0,0;...
         1/tie_a,0,0,0,-1/tie_a,0,0,0;...
         1/tie_n,0,0,0,0,-1/tie_n,0,0;...
         0,1/tei,0,0,0,0,-1/tei,0;...
         0,1/tii,0,0,0,0,0,-1/tii];
    B = [Ie/te;Ii/ti;0;0;0;0;0;0];
    %C = [0,0,((1-q)-qdiff1)*(Jee+fpos)/te,(q+qdiff1)*(Jee+fpos)/te,0,0,-Jei/te,0];
    C = [1,0,0,0,0,0,0,0];
    D = 0;

    sys1 = ss(A,B,C,D);
    S = stepinfo(sys1);
    %step(sys1,t)
    rise = S.RiseTime;
 end


 function sys1 = system(diff,q,w,k)
    te  = 0.02;
    ti  = 0.01;
    tee_a = 0.005;
    tee_n = 0.1;
    tie_a = 0.005;
    tie_n = 0.1;
    tei = 0.01;
    tii = 0.01;
    qdiff1 = diff;
    qdiff2 = 0;
    Jee = w;
    Jei = k*w;
    Jie = w;
    Jii = k*w;
    fpos = 0;

    A = [-1/te,0,((1-q)-qdiff1)*(Jee+fpos)/te,(q+qdiff1)*(Jee+fpos)/te,0,0,-Jei/te,0;...
         0,-1/ti,0,0,((1-q)-qdiff2)*Jie/ti,(q+qdiff2)*Jie/ti,0,-Jii/ti;...
         1/tee_a,0,-1/tee_a,0,0,0,0,0;...
         1/tee_n,0,0,-1/tee_n,0,0,0,0;...
         1/tie_a,0,0,0,-1/tie_a,0,0,0;...
         1/tie_n,0,0,0,0,-1/tie_n,0,0;...
         0,1/tei,0,0,0,0,-1/tei,0;...
         0,1/tii,0,0,0,0,0,-1/tii];
    B = [Ie/te;Ii/ti;0;0;0;0;0;0];
    C = [1,0,0,0,0,0,0,0];
    D = 0;

    sys1 = ss(A,B,C,D);
 end

end


















