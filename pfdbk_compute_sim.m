 function [sol] = pfdbk_compute_sim(w,q,I_step,w_in,q_in,Io_end,t_end)

te  = 0.02;
tee_a = 0.005;
tee_n = 0.1;

tau = [0;0;0;1/tee_a;1/tee_n];


t0 = 0;
y0 = [0;0;0;0;0]; 

A = [-1/te,w*(1-q)/te,w*q/te,(1-q_in)/te,q_in/te;...
    1/tee_a,-1/tee_a,0,0,0;...
    1/tee_n,0,-1/tee_n,0,0;...
    0,0,0,-1/tee_a,0;...
    0,0,0,0,-1/tee_n];

tspan = [t0 t_end]; 
options = odeset('MaxStep',0.001);
sol = ode45(@odefun,tspan,y0,options);



function dy = odefun(t,y)
    
    input = step_input(t);
    
    dy =  A * y + tau .* input ;
end



function input = step_input(t)
    if t<=Io_end
        input = [0;0;0;I_step;I_step];
    else
        input = [0;0;0;0;0];
    end
end




end


















