function [rise_time] = pfdbk_compute_RT(w,q,q_in)

rise_time = risetime(q,q_in,w);

function rise = risetime(q,q_in,w)
    te  = 0.02;
    tee_a = 0.005;
    tee_n = 0.1;

    A = [-1/te,(1-q)*w/te,q*w/te,(1-q_in)/te,q_in/te;...
         1/tee_a,-1/tee_a,0,0,0;...
         1/tee_n,0,-1/tee_n,0,0;...
         0,0,0,-1/tee_a,0;...
         0,0,0,0,-1/tee_n];
    B = [0;0;0;1/tee_a;1/tee_n];
    C = [1,0,0,0,0];
    D = 0;

    sys1 = ss(A,B,C,D);
    S = stepinfo(sys1);
    %step(sys1,t)
    rise = S.RiseTime;
end



end


















