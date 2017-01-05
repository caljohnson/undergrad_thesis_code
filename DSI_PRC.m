function [  ] = DSI_PRC( I, delV, tao, taoS, epsilon )
%Calculates several PRCs for the given parameters
%   Finds the stable period, then gives epsilon input and finds PRC
%   analytically

%steady-state constants
[T , stabilityCheck] = periodLogMap(I, delV, tao, taoS, 1);     %stable period

V = I*(1 - exp(-tao)) + delV*exp((-T + tao)/taoS)*(taoS/(taoS - 1))*(exp(-tao/taoS) - exp(-tao));   %stable voltage at which the delay kicks in

%simulation constants
dt = .001;                                   %time step
delTheta = zeros(floor(T/dt) + 1, 1);       %new period difference, single iterate
delThetaTOT = zeros(floor(T/dt) + 1, 1);    %phase difference, eventual iterate
omega = 0:dt:floor(T/dt)*dt;                %array of stimulus timing values
options = optimset('Display','off');        %represses extraneous output


%firing map functions
dF = (1-I)*exp(T - tao) - delV*exp((T-tao)*(taoS-1)/(taoS));        %derivative of F(Tk+1) at T*
dG = delV/(taoS-1)*exp(-T/taoS)*(exp(-tao+tao/taoS) - 1);           %derivative of G(Tk) at T*


%1st Order PRC
for i = 1:tao/dt+1  %for omega between 0 and tao
    delTheta(i) = epsilon*exp(-tao+omega(i))/dF;
    delThetaTOT(i) = delThetaTOT(i) + delTheta(i);
    %disp(delTheta(i));
    
end
for i = floor(tao/dt + 2) : floor(T/dt)   %for omega between tao and T*
   Vomg = I*(1-exp(-omega(i)+tao)) + V*exp(-omega(i)+tao) + delV*(taoS/(taoS-1))*(exp((-omega(i)+tao)/taoS) - exp(-omega(i)+tao));
   To = fsolve(@(x) (-1 + I*(1-exp(-x+omega(i))) + (Vomg+epsilon)*exp(-x+omega(i)) + delV*(taoS/(taoS-1))*(exp((-x+tao)/taoS) - exp(-x + omega(i) - omega(i)/taoS + tao/taoS))), T, options); %verified right equation
   delTheta(i) = To - T;
   delThetaTOT(i) = delThetaTOT(i) + delTheta(i);
end

 %figure(2); clf; plot(omega/T, -delTheta/epsilon/T, 'r');

%2nd Order PRC
for i = 1:T/dt + 1
    delTheta(i) = (dG/dF)*delTheta(i);
    delThetaTOT(i) = delThetaTOT(i) + delTheta(i);
end

%figure(8); clf; plot(omega, -delTheta/epsilon, 'r');
%3rd Order PRC
for i = 1:T/dt + 1
    delTheta(i) = (dG/dF)*delTheta(i);
    delThetaTOT(i) = delThetaTOT(i) + delTheta(i);
end

figure(9); clf; plot(omega, -delTheta/epsilon, 'r');

%Eventual PRC
for j = 1:5
    for i = 1:T/dt + 1
     delTheta(i) = (dG/dF)*delTheta(i);
     delThetaTOT(i) = delThetaTOT(i) + delTheta(i);
    end
end

figure(9); clf; plot(omega/T, -delThetaTOT/(epsilon*T), 'r');

num_DSI_PRC(T, I, delV, tao, taoS, epsilon);
disp(stabilityCheck);
end
