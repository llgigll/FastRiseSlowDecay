function rise = pfdbk_std_impact_linrisetime(w)

q = 0.5;
rise = zeros(length(w),1);

for i = 1:length(w)
    rise(i) = risetime(q, w(i));
end

function rise = risetime(q,w)
    te  = 0.02;
    tee_a = 0.005;
    tee_n = 0.1;

    A = [-1/te,(1-q)*w/te,q*w/te;...
         1/tee_a,-1/tee_a,0;...
         1/tee_n,0,-1/tee_n];
    B = [1/te;0;0];
    C = [1,0,0];
    D = 0;

    sys1 = ss(A,B,C,D);
    S = stepinfo(sys1);
    %step(sys1,t)
    rise = S.RiseTime;
end


end


















