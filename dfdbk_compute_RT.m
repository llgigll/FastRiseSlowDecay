function [rise_time] = dfdbk_compute_RT(w,q,q_in,diff,k)

rise_time = risetime(diff,q,q_in,w,k);


 function rise = risetime(diff,q,q_in,w,k)
    te  = 0.02;
    ti  = 0.01;
    tee_a = 0.005;
    tee_n = 0.1;
    tie_a = 0.005;
    tie_n = 0.1;
    tei = 0.01;
    tii = 0.01;
    Jee = w;
    Jei = k*w;
    Jie = w;
    Jii = k*w;

    A = [-1/te,0,(1-q-diff)*Jee/(te),(q+diff)*Jee/(te),0,0,-Jei/te,0,(1-q_in)/te,q_in/te;...
         0,-1/ti,0,0,(1-q)*Jie/(ti),q*Jie/(ti),0,-Jii/ti,0,0;...
         1/tee_a,0,-1/tee_a,0,0,0,0,0,0,0;...
         1/tee_n,0,0,-1/tee_n,0,0,0,0,0,0;...
         1/tie_a,0,0,0,-1/tie_a,0,0,0,0,0;...
         1/tie_n,0,0,0,0,-1/tie_n,0,0,0,0;...
         0,1/tei,0,0,0,0,-1/tei,0,0,0;...
         0,1/tii,0,0,0,0,0,-1/tii,0,0;...
         0,0,0,0,0,0,0,0,-1/tee_a,0;...
         0,0,0,0,0,0,0,0,0,-1/tee_n];
    B = [0;0;0;0;0;0;0;0;1/tee_a;1/tee_n];
    C = [1,0,0,0,0,0,0,0,0,0];
    D = 0;

    sys1 = ss(A,B,C,D);
    S = stepinfo(sys1);
    %step(sys1,t)
    rise = S.RiseTime;
 end



end


















