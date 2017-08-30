function [sol,max_eigs,ss_all] = dfdbk_std_split_sim(w,k,q,diff,tr1_e,tr2_e,U1_e,U2_e,I_step,w_in,q_in,Io_end,t_end)
% Simulte the response of the derivative feedback network with STD on the
% excitatory connections.

Io_begin = 0;
Io_step = 0.0;
I = I_step;

t0 = 0;

Jee = w;
Jei = k*w;
Jie = w;
Jii = k*w;

te  = 0.02;
ti  = 0.01;
tee_a = 0.005;
tee_n = 0.1;
tie_a = 0.005;
tie_n = 0.1;
tei = 0.01;
tii = 0.01;
tin_a = 0.005;
tin_n = 0.1;

qie = q;
qee = q;
tau = [zeros(8,1);1/tin_a;1/tin_n;1/tr1_e;1/tr2_e];

max_eigs  = dfdbk_std_split_steadystate(w,k,q,diff,tr1_e,tr2_e,U1_e,U2_e,I_step)


A = [-1/te,0,(1-qee-diff)*Jee/(te),(qee+diff)*Jee/(te),0,0,-Jei/te,0,(1-q_in)*w_in/te,q_in*w_in/te,0,0;...
     0,-1/ti,0,0,(1-qie)*Jie/(ti),qie*Jie/(ti),0,-Jii/ti,0,0,0,0;...
     1/tee_a,0,-1/tee_a,0,0,0,0,0,0,0,0,0;...
     1/tee_n,0,0,-1/tee_n,0,0,0,0,0,0,0,0;...
     1/tie_a,0,0,0,-1/tie_a,0,0,0,0,0,0,0;...
     1/tie_n,0,0,0,0,-1/tie_n,0,0,0,0,0,0;...
     0,1/tei,0,0,0,0,-1/tei,0,0,0,0,0;...
     0,1/tii,0,0,0,0,0,-1/tii,0,0,0,0;...
     0,0,0,0,0,0,0,0,-1/tin_a,0,0,0;...
     0,0,0,0,0,0,0,0,0,-1/tin_n,0,0;...
     0,0,0,0,0,0,0,0,0,0,-1/tr1_e,0;...
     0,0,0,0,0,0,0,0,0,0,0,-1/tr2_e];

 
ss_all =  steady_states(w,k,q,diff,tr1_e,tr2_e,U1_e,U2_e,I_step);
tspan = [t0 t_end];

options = odeset('MaxStep',0.001);
y0 = [zeros(10,1);1;1];
sol = ode45(@odefun,tspan,y0,options);

function dy = odefun(t,y)
    
    input = step_input(t);
    
    B = A;
    
    B(3,1) = A(3,1)*y(11);
    B(4,1) = A(4,1)*y(12);
    B(5,1) = A(5,1)*y(12);
    B(6,1) = A(6,1)*y(11);
    
    B(11,11) = A(11,11)-U1_e*y(1);
    B(12,12) = A(12,12)-U2_e*y(1);
    dy =  B * y + tau .* input ;
    
end



function input = step_input(t)
    if t >= Io_begin && t <= Io_step
        input = [0;0;0;0;0;0;0;0;I_step;I_step;1;1];
    elseif t > Io_step && t<=Io_end
        input = [0;0;0;0;0;0;0;0;I;I;1;1];
    else
        input = [0;0;0;0;0;0;0;0;0;0;1;1];
    end
end

function [ss_all] =  steady_states(w,k,q,diff,tr1_e,tr2_e,U1_e,U2_e,I_step)
    third =   tr1_e*U1_e*tr2_e*U2_e;
    second =  (tr1_e*U1_e+tr2_e*U2_e)-(1-q-diff-(k*w)./(1+k*w)*q)*w*tr2_e*U2_e...
        -(q+diff-(k*w)./(1+k*w)*(1-q))*w*tr1_e*U1_e-tr1_e*U1_e*tr2_e*U2_e*I_step;
    first = 1-(1-q-diff-(k*w)./(1+k*w)*q)*w-(q+diff-(k*w)./(1+k*w)*(1-q))*w...
        -(tr1_e*U1_e+tr2_e*U2_e)*I_step;
    zeroeth = -I_step;
    Re_coeffs = [third,second,first,zeroeth];
    roots(Re_coeffs)
    Re_ss = max(roots(Re_coeffs));
    
    x1_ss = 1/(1+tr1_e*U1_e*Re_ss);
    x2_ss = 1/(1+tr2_e*U2_e*Re_ss);
    
    See_a = x1_ss*Re_ss;
    See_n = x2_ss*Re_ss;
    Sie_a = x2_ss*Re_ss;
    Sie_n = x1_ss*Re_ss;
    
    Ri_ss = ((1-q)*x2_ss+q*x1_ss)*w*Re_ss./(1+k*w);
    
    Sii = Ri_ss;
    Sei = Ri_ss;
    
    Sin_a = I_step;
    Sin_n = I_step;
    
    ss_all = [Re_ss,Ri_ss,See_a,See_n,Sie_a,Sie_n,Sei,Sii,Sin_a,Sin_n,x1_ss,x2_ss];
end



end


















