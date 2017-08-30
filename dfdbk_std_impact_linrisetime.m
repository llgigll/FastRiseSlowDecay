function rise = dfdbk_std_impact_linrisetime(delta_q,w,k)


q = 0.5;
rise = zeros(length(w),1);
for i = 1:length(w)
    rise(i) = risetime(delta_q, q, w(i), k);
end



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
    B = [1/te;0/ti;0;0;0;0;0;0];
    %C = [0,0,((1-q)-qdiff1)*(Jee+fpos)/te,(q+qdiff1)*(Jee+fpos)/te,0,0,-Jei/te,0];
    C = [1,0,0,0,0,0,0,0];
    D = 0;

    sys1 = ss(A,B,C,D);
    S = stepinfo(sys1);
    %step(sys1,t)
    rise = S.RiseTime;
end


end


















