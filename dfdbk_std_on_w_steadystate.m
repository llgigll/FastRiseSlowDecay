function max_eigs  = dfdbk_std_on_w_steadystate(w,k,q,qdiff,tr_e,U_e,Ie)
% What about if the steady state for Ie=0 is not zero??????????

Re_coeffs = [1,(1+w*(k-1))./(tr_e*U_e*(1+k*w))-Ie,-Ie./(tr_e*U_e)];
Re_ss = max(roots(Re_coeffs));
x_ss = 1/(1+tr_e*U_e*Re_ss);
% See_a = x_ss*Re_ss;
% See_n = x_ss*Re_ss;
% Sie_a = x_ss*Re_ss;
% Sie_n = x_ss*Re_ss;
% Ri_ss = x_ss*w*Re_ss./(1+k*w);
% Sii = Ri_ss;
% Sei = Ri_ss;
% Sin_a = I_step;
% Sin_n = I_step;


te  = 0.02;
ti  = 0.01;
tee_a = 0.005;
tee_n = 0.1;
tie_a = 0.005;
tie_n = 0.1;
tei = 0.01;
tii = 0.01;



A_high = [-1/te,0,((1-q)-qdiff)*w/te,(q+qdiff)*w/te,0,0,-k*w/te,0,0;...
    0,-1/ti,0,0,(1-q)*w/ti,q*w/ti,0,-k*w/ti,0;...
    x_ss/tee_a,0,-1/tee_a,0,0,0,0,0,Re_ss/tee_a;...
    x_ss/tee_n,0,0,-1/tee_n,0,0,0,0,Re_ss/tee_n;...
    x_ss/tie_a,0,0,0,-1/tie_a,0,0,0,Re_ss/tie_a;...
    x_ss/tie_n,0,0,0,0,-1/tie_n,0,0,Re_ss/tie_n;...
    0,1/tei,0,0,0,0,-1/tei,0,0;...
    0,1/tii,0,0,0,0,0,-1/tii,0;...
    -U_e*x_ss,0,0,0,0,0,0,0,-1/tr_e-U_e*Re_ss];

A_low = [-1/te,0,((1-q)-qdiff)*w/te,(q+qdiff)*w/te,0,0,-k*w/te,0,0;...
    0,-1/ti,0,0,(1-q)*w/ti,q*w/ti,0,-k*w/ti,0;...
    1/tee_a,0,-1/tee_a,0,0,0,0,0,0;...
    1/tee_n,0,0,-1/tee_n,0,0,0,0,0;...
    1/tie_a,0,0,0,-1/tie_a,0,0,0,0;...
    1/tie_n,0,0,0,0,-1/tie_n,0,0,0;...
    0,1/tei,0,0,0,0,-1/tei,0,0;...
    0,1/tii,0,0,0,0,0,-1/tii,0;...
    -U_e,0,0,0,0,0,0,0,-1/tr_e];

A_lin_high = [-1/te,0,((1-q)-qdiff)*x_ss*w/te,(q+qdiff)*x_ss*w/te,0,0,-k*w/te,0;...
    0,-1/ti,0,0,(1-q)*x_ss*w/ti,q*x_ss*w/ti,0,-k*w/ti;...
    1/tee_a,0,-1/tee_a,0,0,0,0,0;...
    1/tee_n,0,0,-1/tee_n,0,0,0,0;...
    1/tie_a,0,0,0,-1/tie_a,0,0,0;...
    1/tie_n,0,0,0,0,-1/tie_n,0,0;...
    0,1/tei,0,0,0,0,-1/tei,0;...
    0,1/tii,0,0,0,0,0,-1/tii];

A_lin_low = [-1/te,0,((1-q)-qdiff)*w/te,(q+qdiff)*w/te,0,0,-k*w/te,0;...
    0,-1/ti,0,0,(1-q)*w/ti,q*w/ti,0,-k*w/ti;...
    1/tee_a,0,-1/tee_a,0,0,0,0,0;...
    1/tee_n,0,0,-1/tee_n,0,0,0,0;...
    1/tie_a,0,0,0,-1/tie_a,0,0,0;...
    1/tie_n,0,0,0,0,-1/tie_n,0,0;...
    0,1/tei,0,0,0,0,-1/tei,0;...
    0,1/tii,0,0,0,0,0,-1/tii];

eig(A_high)
nonlin_stab_high = max(real(eig(A_high)));
nonlin_stab_low = max(real(eig(A_low)));
lin_stab_high = max(real(eig(A_lin_high)));
lin_stab_low = max(real(eig(A_lin_low)));

max_eigs = [nonlin_stab_high,nonlin_stab_low,lin_stab_high,lin_stab_low];
        
end



























