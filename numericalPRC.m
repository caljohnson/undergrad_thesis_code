%Numerical PRC for a single DSI cell modeled by LIF
clear all;

%Neuron properties- These are actually unused in the following computation
g = 0.2;        % gL, resistance constant for leakage channel  
E = 0.1;        % EL, voltage potential across leakage channel
I = 0.1;        % I, constant input current
Vtho = 0.5;     %Vthreshold, original equation
Vreso = 0.1;    %Vreset, original

%delay properties
delV = -0.5;    %weight of self-inhibition
tao = 1.1;      %delay in seconds
taoS = 2;       %synaptic time constant

%dimensionless system properties
Vthn = 1;                           %scaled threshold
Vresn = 0;                          %scaled reset value
I = 1.1;                            %I' fixed for computation

%start the simulation at the steady state

%steady-state constants
Tstar = 5.5708;     %stable period
Vstar = I*(1 - exp(-tao)) + delV*exp((-Tstar + tao)/taoS)*(taoS/(taoS - 1))*(exp(-tao/taoS) - exp(-tao));   %stable voltage at which the delay kicks in

%simulation constants
dt = 0.0001;         %time step
tF = round(6*Tstar);%final time
Last = tF/dt + 1;   %final time, scaled to unit size
t = 0:dt:tF;

%input stimulus
epsilon = 0.001;       %weight of input stimulus
do = .001;             %omega step            

omega = 0: do: round(Tstar) - do;       %initialize values of omega array
delTheta = zeros(round(Tstar)/do,7);    %initialize PRC array values to 0

%overarching for loop through values of omega
for o = 1:round(Tstar/do)

%begin simulation for omega = omega(o)

%initial conditions
V(1) = Vresn;                   %init cond for DSI cell
S(1) = exp((-Tstar+tao)/taoS);  %start S at steady state, value at Tstar
n = 2;                          %keeps track of pointer in V and S vectors
m = tao;                        %default delay variable to tao seconds
delay = 1;                      %default initial delay trigger to on
PRCorder = 1;                   %phase response curve order, begins with first
Vnew = V(1);                    %start voltage at 0
Snew = S(1);                    %start S at steady state, value at Tstar

%iterate through time steps to observe system behavior
for x=2:1:Last
    Vold = Vnew;                %keeps track of voltage from previous step  
    Sold = Snew;                %keeps track of old synaptic current from previous step
    
    if(t(x) == omega(o))        %if the current timing corresponds to the stimulus timing
       Vold = Vold + epsilon;       %then we change the voltage by amount epsilon
    end
    
    if ( delay == 1 && m <= 0)  %if the delay has been triggered and the tao-second countdown has expired
       Sold = 1;                    %then jump the synaptic current up to 1
       delay = 0;                       %and turn off the delay trigger
    end
    
    if(Vold >= Vthn)        %when the voltage hits threshold
        Vnew = 0;               %reset voltage to 0
        m = tao;                %start tao-second timer
        delay = 1;              %turn on delay trigger
        
                                %calculate PRC of order PRCorder
        delTheta(o, PRCorder) = PRCorder*Tstar - t(x);      %find total phase change for PRCorder
        if(PRCorder>1)                                      %if the order of the PRC is greater than 1
           for i = 1: PRCorder - 1                              %then subtract all lesser order phase shifts
               delTheta(o, PRCorder) = delTheta(o, PRCorder) - delTheta(o, PRCorder-i); 
           end                                                  %to obtain the individaul phase shift of order PRCorder
        end
        
        PRCorder = PRCorder+1;  %up the order of the PRC
           
           
    else                                            %if no other event happens, then simply
        Vnew = Vold + (-Vold + I + delV*Sold)*dt;       %change the voltage according to the DE
        Snew = Sold - (Sold/taoS)*dt;                   %change the synaptic current according to the DE
        if(m>0)                                         %and if the delay timer is active
            m = m - dt;                                     %then count it down by the time step
        end
    end

    
    if(mod(x,10) == 0)  %only store every tenth voltage to plot later
        V(n) = Vnew;
        S(n) = Snew;
        n = n + 1;
    end
end

delTheta(o, 7) = delTheta(o,1) + delTheta(o,2) + delTheta(o,3) + delTheta(o, 4) + delTheta(o, 5) + delTheta(o, 6);  %calculate eventual phase shift

end
 figure(1); clf; PS = plot(V, '-b'); hold on; PR = plot(S, 'r'); %plot voltage and synaptic current
 figure(2); clf; P1 = plot(omega, delTheta(:, 1), '-b');                %plot 1st order PRC   
 figure(3); clf; P2 = plot(omega, delTheta(:, 2), '-b');                %plot 2nd order PRC
 figure(4); clf; P3 = plot(omega, delTheta(:, 3), '-b');                %plot 3rd order PRC
 figure(5); clf; PE = plot(omega, delTheta(:, 7), '-b');                %plot Eventual PRC
