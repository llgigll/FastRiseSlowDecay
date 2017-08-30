function [rise_time, decay_time] = dfdbk_std_manysynapses_RandD(w,k,q,qdiff,tr,U,I_step,w_in,q_in,stim_setup)
% Simulate the response of the derivative feedback network with STD on the
% excitatory connections.

stim = stim_setup(1);
isi = stim_setup(2);
num_reps = stim_setup(3);
t_end = (stim+isi)*num_reps;
I = I_step;


N = length(U); % Number of different synapses on each projection

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



% For used for constant values or external inputs
tau = [0;0;1/tee_a;1/tee_n;1/tie_a;1/tie_n;0;0;1/tin_a;1/tin_n];

A = [-1/te,0,Jee/te,Jee/te,0,0,-Jei/te,0,(1-q_in)*w_in/te,q_in*w_in/te;...
     0,-1/ti,0,0,Jie/ti,Jie/ti,0,-Jii/ti,0,0;...
     0,0,-1/tee_a,0,0,0,0,0,0,0;...
     0,0,0,-1/tee_n,0,0,0,0,0,0;...
     0,0,0,0,-1/tie_a,0,0,0,0,0;...
     0,0,0,0,0,-1/tie_n,0,0,0,0;...
     0,1/tei,0,0,0,0,-1/tei,0,0,0;...
     0,1/tii,0,0,0,0,0,-1/tii,0,0;...
     0,0,0,0,0,0,0,0,-1/tin_a,0;...
     0,0,0,0,0,0,0,0,0,-1/tin_n];

ss_all =  steady_states(w,k,q,qdiff,tr,U,I_step);

% compute rise time of the network
events = @(t,y) event_types(t,y,ss_all(1),true);
tspan = [t0 t_end];
options = odeset('MaxStep',0.001,'nonnegative',1:12,'events',events);
y0 = [zeros(10,1);ones(N,1)];
sol = ode45(@odefun,tspan,y0,options);
ye = sol.y(1,:);
T = sol.x;
rise_time = T(end)-T(find(ye>0.1*ss_all(1),1,'first'));


I_step = 0.0;
I = I_step;
events = @(t,y) event_types(t,y,ss_all(1),false);
options = odeset('MaxStep',0.001,'events',events);
y0 = ss_all;
sol = ode45(@odefun,tspan,y0,options);

ye = sol.y(1,:);
T = sol.x;
decay_time = T(end)-T(find(ye<0.9*ss_all(1),1,'first'));

function dy = odefun(t,y)
      
    B = A;
    
    % Compute time dependent input
    input = step_input(t);

    % Pulling out std values
    y_std = y(11:(11+N-1));
    
    % Sum Inputs to synapses:
    % Excitatory to Excitatory
    input(3) = (1/N)*y(1)*dot((1-q-qdiff),y_std);  %AMPA
    input(4) = (1/N)*y(1)*dot((q+qdiff),y_std);     %NMDA
    % Excitatory to Inhibitory
    input(5) = (1/N)*y(1)*dot((1-q),flip(y_std));   %AMPA
    input(6) = (1/N)*y(1)*dot(q,flip(y_std));       %NMDA
    
    % Evaluate change in STD;
    dy_std = (1-y_std)./tr-U.*y_std.*y(1);
    
    
    dy_net =  B * y(1:10) + tau .* input ;
    
    dy = [dy_net;dy_std];
    
end



function input = step_input(t)
    rep_time = mod(t,(isi+stim));
    if rep_time <= stim
        input = [0;0;0;0;0;0;0;0;I;I];
    else
        input = [0;0;0;0;0;0;0;0;0;0];
    end
end


function [ss_all] =  steady_states(w,k,q,qdiff,tr,u,I_step)
    
    f = @(u1,tr1,u2,tr2,w,k,I) roots([1,-I+(1+0.5*(k*w^2)/(1+k*w)-0.5*w)*...
    (1/(u1*tr1)+1/(u2*tr2)),-I*(1/(u1*tr1)+1/(u2*tr2))+...
    (1+(k*w^2)/(1+k*w)-w)/(u1*tr1*u2*tr2),-I/(u1*tr1*u2*tr2)]);

    Re_ss = max(f(u(1),tr(1),u(2),tr(2),w,k,I));
    
    x1_ss = 1/(1+tr(1)*u(1)*Re_ss);
    x2_ss = 1/(1+tr(2)*u(2)*Re_ss);
    See_a = 0.5*((1-q(1)-qdiff)*x1_ss+(1-q(2)-qdiff)*x2_ss)*Re_ss;
    See_n = 0.5*((q(1)+qdiff)*x1_ss+(q(2)+qdiff)*x2_ss)*Re_ss;
    Sie_a = 0.5*((1-q(1))*x2_ss+(1-q(2))*x1_ss)*Re_ss;
    Sie_n = 0.5*(q(1)*x2_ss+q(2)*x1_ss)*Re_ss;
    Ri_ss = 0.5*(x1_ss+x2_ss)*w*Re_ss./(1+k*w);
    Sii = Ri_ss;
    Sei = Ri_ss;
    Sin_a = I_step;
    Sin_n = I_step;
    
    ss_all = [Re_ss,Ri_ss,See_a,See_n,Sie_a,Sie_n,Sei,Sii,Sin_a,Sin_n,x1_ss,x2_ss];
end


function [value,isterminal,direction] = event_types(t,y,Re_ss,rise)

    if rise
        rise_val = y(1)-Re_ss*0.9;
        value = rise_val;
        isterminal = 1;   % Stop at local minimum
        direction = 1;  % [local minimum, local maximum]
    else
        decay_val = y(1)-Re_ss*0.1;
        value = decay_val;
        isterminal = 1;   % Stop at local minimum
        direction = -1;  % [local minimum, local maximum]
    end

end 


end


















