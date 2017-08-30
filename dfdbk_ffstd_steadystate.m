function Re_ss = dfdbk_ffstd_steadystate(w,k,tr,U,tr_ff,U_ff,I_step,w_in)

x_in = 1/(1+tr_ff*U_ff*I_step);
I_in = w_in*I_step*x_in;

f = @(u1,tr1,u2,tr2,w,k,I_in) roots([1,-I_in+(1+0.5*(k*w^2)/(1+k*w)-0.5*w)*...
(1/(u1*tr1)+1/(u2*tr2)),-I_in*(1/(u1*tr1)+1/(u2*tr2))+...
(1+(k*w^2)/(1+k*w)-w)/(u1*tr1*u2*tr2),-I_in/(u1*tr1*u2*tr2)]);

Re_ss = max(f(U(1),tr(1),U(2),tr(2),w,k,I_in));

end


















