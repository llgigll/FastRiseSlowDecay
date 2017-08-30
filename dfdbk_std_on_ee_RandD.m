function [rise_time,decay_time,max_eigs] = dfdbk_std_on_ee_RandD(w,k,q,diff,tr_e,U_e,I_step,w_in,q_in)
% Compute the rise and decay time of the derivative feedback network with
% STD on the excitatory connections.



Io_begin = 0;
Io_step = 0.0;
Io_end = 100.0;

t0 = 0;
t_end = 100.0;

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

qie = q
qee = q
tau = [zeros(8,1);1/tin_a;1/tin_n;1/tr_e];

max_eigs  = dfdbk_std_on_w_steadystate(w,k,q,diff,tr_e,U_e,I_step);


A = [-1/te,0,(1-qee-diff)*Jee/(te),(qee+diff)*Jee/(te),0,0,-Jei/te,0,(1-q_in)*w_in/te,q_in*w_in/te,0;...
     0,-1/ti,0,0,(1-qie)*Jie/(ti),qie*Jie/(ti),0,-Jii/ti,0,0,0;...
     1/tee_a,0,-1/tee_a,0,0,0,0,0,0,0,0;...
     1/tee_n,0,0,-1/tee_n,0,0,0,0,0,0,0;...
     1/tie_a,0,0,0,-1/tie_a,0,0,0,0,0,0;...
     1/tie_n,0,0,0,0,-1/tie_n,0,0,0,0,0;...
     0,1/tei,0,0,0,0,-1/tei,0,0,0,0;...
     0,1/tii,0,0,0,0,0,-1/tii,0,0,0;...
     0,0,0,0,0,0,0,0,-1/tin_a,0,0;...
     0,0,0,0,0,0,0,0,0,-1/tin_n,0;...
     0,0,0,0,0,0,0,0,0,0,-1/tr_e];
 

tspan = [t0 t_end];
I = I_step;
ss_all =  steady_states(w,k,tr_e,U_e,w_in,I_step);

events = @(t,y) event_types(t,y,ss_all(1),true);
options = odeset('MaxStep',0.001,'events',events);
y0 = [zeros(10,1);1];
sol = ode45(@odefun,tspan,y0,options);

ye = sol.y(1,:);
T = sol.x;

% figure('color','w')
% plot(T,ye)
% hold on
% plot([T(1),T(end)],[ss_all(1),ss_all(1)],'k--')
% xlabel('Time (s)')
% ylabel('R_e (Hz)')

rise_time = T(end)-T(find(ye>0.1*ss_all(1),1,'first'));

I_step = 0.0;
I = I_step;
events = @(t,y) event_types(t,y,ss_all(1),false);
options = odeset('MaxStep',0.001,'events',events);
y0 = ss_all;
sol = ode45(@odefun,tspan,y0,options);

ye = sol.y(1,:);
T = sol.x;

% figure('color','w')
% plot(T,ye)
% hold on
% plot([T(1),T(end)],[ss_all(1),ss_all(1)],'k--')
% xlabel('Time (s)')
% ylabel('R_e (Hz)')

decay_time = T(end)-T(find(ye<0.9*ss_all(1),1,'first'));



% ye = sol.y(1,:);
% T = sol.x;
% yi = sol.y(2,:);
% stp = sol.y(11,:);
% Sampa = sol.y(9,:);
% Snmda = sol.y(10,:);
% 
% figure('color','w')
% subplot(2,1,1)
% plot(T,ye)
% hold on
% plot([T(1),T(end)],[ss_all(1),ss_all(1)],'k--')
% xlabel('Time (s)')
% ylabel('R_E (Hz)')
% subplot(2,1,2)
% plot(T,yi)
% hold on
% plot([T(1),T(end)],[ss_all(2),ss_all(2)],'k--')
% xlabel('Time (s)')
% ylabel('R_I (Hz)')


function dy = odefun(t,y)
    
    input = step_input(t);
    
    B = A;
%     B(1,3) = A(1,3)*y(11);
%     B(1,4) = A(1,4)*y(11);
%     B(2,5) = A(2,5)*y(11);
%     B(2,6) = A(2,6)*y(11);
%     B(1,7) = A(1,7)*y(12);
%     B(2,8) = A(2,8)*y(12);
    
    B(3,1) = A(3,1)*y(11);
    B(4,1) = A(4,1)*y(11);
    B(5,1) = A(5,1)*y(11);
    B(6,1) = A(6,1)*y(11);
    
    B(11,11) = A(11,11)-U_e*y(1);
    dy =  B * y + tau .* input ;
    
end



function input = step_input(t)
    if t >= Io_begin && t <= Io_step
        input = [0;0;0;0;0;0;0;0;I_step;I_step;1];
    elseif t > Io_step && t<=Io_end
        input = [0;0;0;0;0;0;0;0;I;I;1];
    else
        input = [0;0;0;0;0;0;0;0;0;0;1];
    end
end

function [ss_all] =  steady_states(w,k,tr_e,U_e,w_in,I_step)
    Re_coeffs = [1,(1+w*(k-1))./(tr_e*U_e*(1+k*w))-w_in*I_step,-w_in*I_step./(tr_e*U_e)];
    Re_ss = max(roots(Re_coeffs));
    x_ss = 1/(1+tr_e*U_e*Re_ss);
    See_a = x_ss*Re_ss;
    See_n = x_ss*Re_ss;
    Sie_a = x_ss*Re_ss;
    Sie_n = x_ss*Re_ss;
    Ri_ss = x_ss*w*Re_ss./(1+k*w);
    Sii = Ri_ss;
    Sei = Ri_ss;
    Sin_a = w_in*I_step;
    Sin_n = w_in*I_step;
    
    ss_all = [Re_ss,Ri_ss,See_a,See_n,Sie_a,Sie_n,Sei,Sii,Sin_a,Sin_n,x_ss];
end

function [value,isterminal,direction] = event_types(t,y,Re_ss,rise)

    if rise
        other_val = 2.0-t;
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


















