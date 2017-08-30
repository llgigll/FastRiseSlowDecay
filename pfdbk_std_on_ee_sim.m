function [sol,ss_all] = pfdbk_std_on_ee_sim(w,q,tr_e,U_e,I_step,w_in,q_in,Io_end,t_end)
% Simulte the response of the derivative feedback network with STD on the
% excitatory connections.

Io_begin = 0;
Io_step = 0.0;

I = I_step;

t0 = 0;

te  = 0.02;
tee_a = 0.005;
tee_n = 0.1;
tin_a = 0.005;
tin_n = 0.1;

tau = [0;0;0;1/tee_a;1/tee_n;1/tr_e];


A = [-1/te,w*(1-q)/te,w*q/te,(1-q_in)*w_in/te,q_in*w_in/te,0;...
    1/tee_a,-1/tee_a,0,0,0,0;...
    1/tee_n,0,-1/tee_n,0,0,0;...
    0,0,0,-1/tin_a,0,0;...
    0,0,0,0,-1/tin_n,0;...
    0,0,0,0,0,-1/tr_e];

 
ss_all =  steady_states(w,tr_e,U_e,I_step);
tspan = [t0 t_end];

options = odeset('MaxStep',0.001);
y0 = [zeros(5,1);1];
sol = ode45(@odefun,tspan,y0,options);



function dy = odefun(t,y)
    
    input = step_input(t);
    B = A;
    B(2,1) = A(2,1)*y(6);
    B(3,1) = A(3,1)*y(6);
    
    
    B(6,6)=A(6,6)-U_e*y(1);
    dy =  B * y + tau .* input ;
end


function [ss_all] =  steady_states(w,tr_e,U_e,I_step)
    Re_coeffs = [1,(1-w)./(tr_e*U_e)-I_step,-I_step./(tr_e*U_e)];
    Re_ss = max(roots(Re_coeffs));
    x_ss = 1/(1+tr_e*U_e*Re_ss);
    See_a = x_ss*Re_ss;
    See_n = x_ss*Re_ss;
    Sin_a = I_step;
    Sin_n = I_step;
    
    ss_all = [Re_ss,See_a,See_n,Sin_a,Sin_n,x_ss];
end


function input = step_input(t)
    if t >= Io_begin && t <= Io_step
        input = [0;0;0;I_step;I_step;1];
    elseif t > Io_step && t<=Io_end
        input = [0;0;0;I;I;1];
    else
        input = [0;0;0;0;0;1];
    end
end


end


















