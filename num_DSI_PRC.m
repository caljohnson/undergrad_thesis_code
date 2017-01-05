function [ ] = num_DSI_PRC( T, I, delV, tao, taoS, epsilon)
%Numerical PRC for DSI cell
%   based on given parameters

%start the simulation at the steady state

%steady-state constants
Tstar = T;     %stable period
Vstar = I*(1 - exp(-tao)) + delV*exp((-Tstar + tao)/taoS)*(taoS/(taoS - 1))*(exp(-tao/taoS) - exp(-tao));   %stable voltage at which the delay kicks in
Vthn = 1;                           %scaled threshold
Vresn = 0;                          %scaled reset value

%simulation constants
dt = 0.001;         %time step
tF = round(6*Tstar);%final time
Last = tF/dt + 1;   %final time, scaled to unit size
% t = zeros(tF/dt+1);

%input stimulus
do = .01;             %omega step            

omega = 0: do: round(Tstar) - do;       %initialize values of omega array
delTheta = zeros(round(Tstar)/do,7);    %initialize PRC array values to 0
Tstars = zeros(7);

%overarching for loop through values of omega


%begin simulation for omega = omega(o)

%initial conditions
V(1) = Vresn;                   %init cond for DSI cell
S(1) = exp((-Tstar+tao)/taoS);  %start S at steady state, value at Tstar
n = 2;                          %keeps track of pointer in V and S vectors
m = tao;                        %default delay variable to tao seconds
delay = 1;                      %default initial delay trigger to on
PRCorder = 1;                   %phase response curve order, begins with first
Vold = V(1);                    %start voltage at 0
Sold = S(1);                    %start S at steady state, value at Tstar
t(1)=0;
%iterate through time steps to observe system behavior
for x=2:1:Last
    
    t(x)=t(x-1)+dt;
    Vnew = Vold + (-Vold + I + delV*Sold)*dt;       %change the voltage according to the DE
    Snew = Sold - (Sold/taoS)*dt;                   %change the synaptic current according to the DE
    if(m>0)                                         %and if the delay timer is active
      m = m - dt;                                     %then count it down by the time step
    end   
    
    if ( delay==1 && m<= 0)  %if the delay has been triggered and the tao-second countdown has expired
       Snew = 1;                    %then jump the synaptic current up to 1
       delay = 0;                       %and turn off the delay trigger
    end
  
    
    if(Vnew >= Vthn)        %when the voltage hits threshold
        m = tao;                %start tao-second timer
        delay = 1;              %turn on delay trigger
        
                                %calculate PRC of order PRCorder
        dtstar = dt*(Vthn-Vold)/(Vnew-Vold);  
        t(x)=t(x-1)+dtstar;
        Tstars(PRCorder) = t(x);      %find precise firing time
        PRCorder = PRCorder+1;  %up the order of the PRC
        Vnew = 0;               %reset voltage to 0
        
        Snew = Sold - (Sold/taoS)*dtstar;                   %change the synaptic current according to the DE
           
    end

    
    if(mod(x,10) == 0)  %only store every tenth voltage to plot later
        V(n) = Vnew;
        S(n) = Snew;
        n = n + 1;
    end
    
    Vold = Vnew;                %keeps track of voltage from previous step  
    Sold = Snew;                %keeps track of old synaptic current from previous step
    
end
disp(PRCorder);
disp(Tstars(1));
%figure(3); clf; plot(V, '-b'); 
%hold on; plot(S, 'r'); 
axis([-Inf, Inf, -Inf, Inf]); %plot voltage and synaptic current


for o = 1:round(Tstar/do)
%initial conditions
V(1) = Vresn;                   %init cond for DSI cell
S(1) = exp((-Tstar+tao)/taoS);  %start S at steady state, value at Tstar
n = 2;                          %keeps track of pointer in V and S vectors
m = tao;                        %default delay variable to tao seconds
delay = 1;                      %default initial delay trigger to on
PRCorder = 1;                   %phase response curve order, begins with first
Vold = V(1);                    %start voltage at 0
Sold = S(1);                    %start S at steady state, value at Tstar
t(1)=0;
for x=2:1:Last
                                        
        t(x)=t(x-1)+dt;
        Vnew = Vold + (-Vold + I + delV*Sold)*dt;       %change the voltage according to the DE
        Snew = Sold - (Sold/taoS)*dt;                   %change the synaptic current according to the DE
        if(m>0)                                         %and if the delay timer is active
            m = m - dt;                                     %then count it down by the time step
        end
    
    if(t(x) <= (omega(o) + dt) && t(x) >= omega(o) )        %if the current timing corresponds to the stimulus timing
       Vnew = Vnew + epsilon;       %then we change the voltage by amount epsilon
    end
    
    if ( delay==1 && m<= 0)  %if the delay has been triggered and the tao-second countdown has expired
       Snew = 1;                    %then jump the synaptic current up to 1
       delay = 0;                       %and turn off the delay trigger
    end
    
    if(Vnew >= Vthn)        %when the voltage hits threshold
        
        m = tao;                %start tao-second timer
        delay = 1;              %turn on delay trigger
        
        dtstar = dt*(Vthn-Vold)/(Vnew-Vold);  
        t(x)=t(x-1)+dtstar;
                                %calculate PRC of order PRCorder
        delTheta(o, PRCorder) = Tstars(PRCorder) - t(x);      %find total phase change for PRCorder
        
        PRCorder = PRCorder+1;  %up the order of the PRC
        Vnew = 0;               %reset voltage to 0
        Snew = Sold - (Sold/taoS)*dtstar;                   %change the synaptic current according to the DE
           
           

    end

    
    if(mod(x,100) == 0)  %only store every tenth voltage to plot later
        V(n) = Vnew;
        S(n) = Snew;
        n = n + 1;
    end
    Vold = Vnew;                %keeps track of voltage from previous step  
    Sold = Snew;                %keeps track of old synaptic current from previous step
end
end


 %figure(2); clf; hold on; plot(omega/T, delTheta(:, 1)/epsilon/T, '-b'); hold on;               %plot 1st order PRC   
 % figure(8); hold on; plot(omega, delTheta(:, 2)/epsilon, '-b');                %plot 2nd order PRC
%  figure(7); hold on; plot(omega, delTheta(:, 3)/epsilon, '-b');                %plot 3rd order PRC
 figure(9); hold on; plot(omega/T, delTheta(:, 5)/(epsilon*T), '-y'); hold on;               %plot Asymptotic PRC
end

