function [sol] = dfdbk_compute_sim(w,k,q,diff,I_step,w_in,q_in,Io_end,t_end)
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
tau = [zeros(8,1);1/tin_a;1/tin_n];


A = [-1/te,0,(1-qee-diff)*Jee/(te),(qee+diff)*Jee/(te),0,0,-Jei/te,0,(1-q_in)*w_in/te,q_in*w_in/te;...
     0,-1/ti,0,0,(1-qie)*Jie/(ti),qie*Jie/(ti),0,-Jii/ti,0,0;...
     1/tee_a,0,-1/tee_a,0,0,0,0,0,0,0;...
     1/tee_n,0,0,-1/tee_n,0,0,0,0,0,0;...
     1/tie_a,0,0,0,-1/tie_a,0,0,0,0,0;...
     1/tie_n,0,0,0,0,-1/tie_n,0,0,0,0;...
     0,1/tei,0,0,0,0,-1/tei,0,0,0;...
     0,1/tii,0,0,0,0,0,-1/tii,0,0;...
     0,0,0,0,0,0,0,0,-1/tin_a,0;...
     0,0,0,0,0,0,0,0,0,-1/tin_n];
 

tspan = [t0 t_end];

options = odeset('MaxStep',0.001);
y0 = zeros(10,1);
sol = ode45(@odefun,tspan,y0,options);



function dy = odefun(t,y)
    
    input = step_input(t);
    
    dy =  A * y + tau .* input ;
    
end



function input = step_input(t)
    if t >= Io_begin && t <= Io_step
        input = [0;0;0;0;0;0;0;0;I_step;I_step];
    elseif t > Io_step && t<=Io_end
        input = [0;0;0;0;0;0;0;0;I;I];
    else
        input = [0;0;0;0;0;0;0;0;0;0];
    end
end




end


















