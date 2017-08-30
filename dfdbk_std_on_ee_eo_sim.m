function [sol,ss_all] = dfdbk_std_on_ee_eo_sim(w,k,q,diff,tr_e,U_e,tr_ff,U_ff,I_step,w_in,q_in,Io_end,t_end)
% Simulate the response of the derivative feedback network with STD on the
% excitatory connections and on the feedforward connections.

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
tau = [zeros(8,1);1/tin_a;1/tin_n;1/tr_e;1/tr_ff];

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
     0,0,0,0,0,0,0,0,0,0,-1/tr_e,0;...
     0,0,0,0,0,0,0,0,0,0,0,-1/tr_ff];

 
ss_all =  steady_states(w,k,tr_e,U_e,tr_ff,U_ff,w_in,I_step);
tspan = [t0 t_end];

options = odeset('MaxStep',0.001);
y0 = [zeros(10,1);1;1];
sol = ode45(@odefun,tspan,y0,options);

function dy = odefun(t,y)
    
    input = step_input(t);
    
    B = A;   
    B(3,1) = A(3,1)*y(11);
    B(4,1) = A(4,1)*y(11);
    B(5,1) = A(5,1)*y(11);
    B(6,1) = A(6,1)*y(11);
    
    B(11,11) = A(11,11)-U_e*y(1);
    B(12,12) = A(12,12)-U_ff*input(end-2);
    input(9:10) = input(9:10)*y(12);
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

function [ss_all] =  steady_states(w,k,tr_e,U_e,tr_ff,U_ff,w_in,I_step)
    x_in = 1/(1+tr_ff*U_ff*I_step);
    I_in = w_in*I_step*x_in;
    Re_coeffs = [1,(1+w*(k-1))./(tr_e*U_e*(1+k*w))-I_in,-I_in./(tr_e*U_e)];
    Re_ss = max(roots(Re_coeffs));
    x_ss = 1/(1+tr_e*U_e*Re_ss);
    See_a = x_ss*Re_ss;
    See_n = x_ss*Re_ss;
    Sie_a = x_ss*Re_ss;
    Sie_n = x_ss*Re_ss;
    Ri_ss = x_ss*w*Re_ss./(1+k*w);
    Sii = Ri_ss;
    Sei = Ri_ss;
    Sin_a = x_in*I_step;
    Sin_n = x_in*I_step;
    
    ss_all = [Re_ss,Ri_ss,See_a,See_n,Sie_a,Sie_n,Sei,Sii,Sin_a,Sin_n,x_ss,x_in];
end



end


















