%Figure 9 - PRC's
    %%plots several asymptotic infinitessimal prc's w/ varying params

    
%%%%%%% PRC 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1st iPRC parameters    
I = 1.1;
alpha = .5;
tau = 1.5;
tau_s = 2;

%yields fixed period length T
T = findT_forFigs(I, alpha, tau, tau_s);

%firing map functions
delV = -alpha;
dF = (1-I)*exp(T - tau) - delV*exp((T-tau)*(tau_s-1)/(tau_s));        %derivative of F(Tk+1) at T*
dG = delV/(tau_s-1)*exp(-T/tau_s)*(exp(-tau+tau/tau_s) - 1);           %derivative of G(Tk) at T*
V = I*(1 - exp(-tau)) + delV*exp((-T + tau)/tau_s)*(tau_s/(tau_s - 1))*(exp(-tau/tau_s) - exp(-tau));   %stable voltage at which the delay kicks in

%simulation constants
dt = .01;                                  %time step
delTheta = zeros(floor(T/dt) +1, 1);       %new period difference, single iterate
delThetaTOT = zeros(floor(T/dt) +1, 1);    %phase difference, eventual iterate
omega = 0:dt:floor(T/dt)*dt;                %array of stimulus timing values
options = optimset('Display','off');        %represses extraneous output
epsilon = .001;

%1st order phase change
for i = 1:floor(tau/dt+1)  %for omega between 0 and tau
    delTheta(i) = epsilon*exp(-tau+omega(i))/dF;
    To = delTheta(i) + T;
    delThetaTOT(i) = delThetaTOT(i) + delTheta(i);

end
for i = floor(tau/dt + 2) : floor(T/dt)   %for omega between tau and T*
   Vomg = I*(1-exp(-omega(i)+tau)) + V*exp(-omega(i)+tau) + delV*(tau_s/(tau_s-1))*(exp((-omega(i)+tau)/tau_s) - exp(-omega(i)+tau));
   To = fsolve(@(x) (-1 + I*(1-exp(-x+omega(i))) + (Vomg+epsilon)*exp(-x+omega(i)) + delV*(tau_s/(tau_s-1))*(exp((-x+tau)/tau_s) - exp(-x + omega(i) - omega(i)/tau_s + tau/tau_s))), To+.001, options);
   delTheta(i) = To - T;
   delThetaTOT(i) = delThetaTOT(i) + delTheta(i);
end
%figure(8);  clf; plot(omega/T, -delTheta/(epsilon*T), 'r');

%asymptotic iPRC
for j = 1:5
    for i = 1:T/dt+1
     delTheta(i) = (dG/dF)*delTheta(i);
     delThetaTOT(i) = delThetaTOT(i) + delTheta(i);
    end
end

figure(9); clf; plot(omega/T, -delThetaTOT/(epsilon*T), 'k'); hold on;
%num_DSI_PRC(T, I, delV, tau, tau_s, epsilon);

%%%%%%% PRC 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%2nd iPRC parameters    
I = 1.1;
alpha = .5;
tau = 2;
tau_s = 2;

%yields fixed period length T
T = findT_forFigs(I, alpha, tau, tau_s);

%firing map functions
delV = -alpha;
dF = (1-I)*exp(T - tau) - delV*exp((T-tau)*(tau_s-1)/(tau_s));        %derivative of F(Tk+1) at T*
dG = delV/(tau_s-1)*exp(-T/tau_s)*(exp(-tau+tau/tau_s) - 1);           %derivative of G(Tk) at T*
V = I*(1 - exp(-tau)) + delV*exp((-T + tau)/tau_s)*(tau_s/(tau_s - 1))*(exp(-tau/tau_s) - exp(-tau));   %stable voltage at which the delay kicks in

%simulation constants
dt = .01;                                  %time step
delTheta = zeros(floor(T/dt) +1, 1);       %new period difference, single iterate
delThetaTOT = zeros(floor(T/dt) +1, 1);    %phase difference, eventual iterate
omega = 0:dt:floor(T/dt)*dt;                %array of stimulus timing values
options = optimset('Display','off');        %represses extraneous output
epsilon = .001;

%1st order phase change
for i = 1:tau/dt+1  %for omega between 0 and tau
    delTheta(i) = epsilon*exp(-tau+omega(i))/dF;
    To = delTheta(i) + T;
    delThetaTOT(i) = delThetaTOT(i) + delTheta(i);

end
for i = floor(tau/dt + 2) : floor(T/dt)   %for omega between tau and T*
   Vomg = I*(1-exp(-omega(i)+tau)) + V*exp(-omega(i)+tau) + delV*(tau_s/(tau_s-1))*(exp((-omega(i)+tau)/tau_s) - exp(-omega(i)+tau));
   To = fsolve(@(x) (-1 + I*(1-exp(-x+omega(i))) + (Vomg+epsilon)*exp(-x+omega(i)) + delV*(tau_s/(tau_s-1))*(exp((-x+tau)/tau_s) - exp(-x + omega(i) - omega(i)/tau_s + tau/tau_s))),  To+.001, options);
   delTheta(i) = To - T;
   delThetaTOT(i) = delThetaTOT(i) + delTheta(i);
end
%figure(8);  clf; plot(omega/T, -delTheta/(epsilon*T), 'r');

%asymptotic iPRC
for j = 1:5
    for i = 1:T/dt+1
     delTheta(i) = (dG/dF)*delTheta(i);
     delThetaTOT(i) = delThetaTOT(i) + delTheta(i);
    end
end

figure(9); plot(omega/T, -delThetaTOT/(epsilon*T), 'b'); hold on;
%num_DSI_PRC(T, I, delV, tau, tau_s, epsilon);

%%%%%% PRC 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%3rd iPRC parameters    
I = 1.1;
alpha = .5;
tau = 1.5;
tau_s = 3;

%yields fixed period length T
T = findT_forFigs(I, alpha, tau, tau_s);

%firing map functions
delV = -alpha;
dF = (1-I)*exp(T - tau) - delV*exp((T-tau)*(tau_s-1)/(tau_s));        %derivative of F(Tk+1) at T*
dG = delV/(tau_s-1)*exp(-T/tau_s)*(exp(-tau+tau/tau_s) - 1);           %derivative of G(Tk) at T*
V = I*(1 - exp(-tau)) + delV*exp((-T + tau)/tau_s)*(tau_s/(tau_s - 1))*(exp(-tau/tau_s) - exp(-tau));   %stable voltage at which the delay kicks in

%simulation constants
dt = .01;                                  %time step
delTheta = zeros(floor(T/dt) +1, 1);       %new period difference, single iterate
delThetaTOT = zeros(floor(T/dt) +1, 1);    %phase difference, eventual iterate
omega = 0:dt:floor(T/dt)*dt;                %array of stimulus timing values
options = optimset('Display','off');        %represses extraneous output
epsilon = .001;

%1st order phase change
for i = 1:tau/dt+1  %for omega between 0 and tau
    delTheta(i) = epsilon*exp(-tau+omega(i))/dF;
    To = delTheta(i) + T;
    delThetaTOT(i) = delThetaTOT(i) + delTheta(i);

end
for i = floor(tau/dt + 2) : floor(T/dt)   %for omega between tau and T*
   Vomg = I*(1-exp(-omega(i)+tau)) + V*exp(-omega(i)+tau) + delV*(tau_s/(tau_s-1))*(exp((-omega(i)+tau)/tau_s) - exp(-omega(i)+tau));
   To = fsolve(@(x) (-1 + I*(1-exp(-x+omega(i))) + (Vomg+epsilon)*exp(-x+omega(i)) + delV*(tau_s/(tau_s-1))*(exp((-x+tau)/tau_s) - exp(-x + omega(i) - omega(i)/tau_s + tau/tau_s))),  To+.001, options);
   delTheta(i) = To - T;
   delThetaTOT(i) = delThetaTOT(i) + delTheta(i);
end
%figure(8);  clf; plot(omega/T, -delTheta/(epsilon*T), 'r');

%asymptotic iPRC
for j = 1:5
    for i = 1:T/dt+1
     delTheta(i) = (dG/dF)*delTheta(i);
     delThetaTOT(i) = delThetaTOT(i) + delTheta(i);
    end
end

figure(9); plot(omega/T, -delThetaTOT/(epsilon*T), 'g'); hold on;
%num_DSI_PRC(T, I, delV, tau, tau_s, epsilon);

%%%%%%% PRC 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%4th iPRC parameters    
I = 1.1;
alpha = .7;
tau = 1.5;
tau_s = 2;

%yields fixed period length T
T = findT_forFigs(I, alpha, tau, tau_s);

%firing map functions
delV = -alpha;
dF = (1-I)*exp(T - tau) - delV*exp((T-tau)*(tau_s-1)/(tau_s));        %derivative of F(Tk+1) at T*
dG = delV/(tau_s-1)*exp(-T/tau_s)*(exp(-tau+tau/tau_s) - 1);           %derivative of G(Tk) at T*
V = I*(1 - exp(-tau)) + delV*exp((-T + tau)/tau_s)*(tau_s/(tau_s - 1))*(exp(-tau/tau_s) - exp(-tau));   %stable voltage at which the delay kicks in

%simulation constants
dt = .01;                                  %time step
delTheta = zeros(floor(T/dt) +1, 1);       %new period difference, single iterate
delThetaTOT = zeros(floor(T/dt) +1, 1);    %phase difference, eventual iterate
omega = 0:dt:floor(T/dt)*dt;                %array of stimulus timing values
options = optimset('Display','off');        %represses extraneous output
epsilon = .001;

%1st order phase change
for i = 1:tau/dt+1  %for omega between 0 and tau
    delTheta(i) = epsilon*exp(-tau+omega(i))/dF;
    To = delTheta(i) + T;
    delThetaTOT(i) = delThetaTOT(i) + delTheta(i);

end
for i = floor(tau/dt + 2) : floor(T/dt)   %for omega between tau and T*
   Vomg = I*(1-exp(-omega(i)+tau)) + V*exp(-omega(i)+tau) + delV*(tau_s/(tau_s-1))*(exp((-omega(i)+tau)/tau_s) - exp(-omega(i)+tau));
   To = fsolve(@(x) (-1 + I*(1-exp(-x+omega(i))) + (Vomg+epsilon)*exp(-x+omega(i)) + delV*(tau_s/(tau_s-1))*(exp((-x+tau)/tau_s) - exp(-x + omega(i) - omega(i)/tau_s + tau/tau_s))),  To+.001, options);
   delTheta(i) = To - T;
   delThetaTOT(i) = delThetaTOT(i) + delTheta(i);
end
%figure(8);  clf; plot(omega/T, -delTheta/(epsilon*T), 'r');

%asymptotic iPRC
for j = 1:5
    for i = 1:T/dt+1
     delTheta(i) = (dG/dF)*delTheta(i);
     delThetaTOT(i) = delThetaTOT(i) + delTheta(i);
    end
end

figure(9); plot(omega/T, -delThetaTOT/(epsilon*T), 'r'); hold on;
%legend('iPRC 1', 'iPRC 2', 'iPRC 3', 'iPRC 4', 'Location', 'northwest');
xlabel('\Omega'); ylabel('\Delta \Theta'); title('Asymptotic iPRCs for DSI cells of varying parameters');

%num_DSI_PRC(T, I, delV, tau, tau_s, epsilon);