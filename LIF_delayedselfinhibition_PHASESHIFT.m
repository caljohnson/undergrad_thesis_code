%A LIF Model for a Single Cell with Delayed Self-Inhibition using
%Exponentially-Decaying Synapses


%This should show us a clear picture of a phase shift caused by a stimulus
%at phase omega
clear; 

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
J = (E - Vreso + I/g)/(Vtho-Vreso); %I' really, dimensionless constant
J = 1.1;                            %I' fixed for computation
T = log(J/(J-1));                   %original period

%simulation constants
dt = 0.001;                 %time step
tF = round(20*T);           %final time
Last = tF/dt + 1;           %final time, scaled to unit size
t = 0:dt:dt*(Last - 1);

%begin simulation

%initial conditions
V(1) = Vresn;              %init cond for DSI cell 
n = 2;                     %keeps track of pointer in V and S vectors
m = tao;                   %default delay variable to tao seconds
delay = 1;                 %default initial delay trigger to on
Vnew = V(1);               %start voltage at 0
Snew = 0;                  %start synaptic current at 0

%iterate through time steps to observe system behavior
for x=2:1:Last
    Vold = Vnew;           %keeps track of voltage from previous step  
    Sold = Snew;           %keeps track of old synaptic current from previous step

    if ( delay == 1 && m <= 0)  %if the delay has been triggered and the tao-second countdown has expired
       Sold = 1;                    %then jump the synaptic current up to 1
       delay = 0;                   %and turn off the delay trigger
    end
    if(Vold >= Vthn)        %when the voltage hits threshold
        Vnew = 0;               %reset voltage to 0
        m = tao;                %start tao-second timer
        delay = 1;              %turn on delay trigger
    
    else                                            %if no other event happens, then simply
        Vnew = Vold + (-Vold + J + delV*Sold)*dt;       %change the voltage according to the DE
        Snew = Sold - (Sold/taoS)*dt;                   %change the synaptic current according to the exponential decay DE
        if(m>0)                                         %and if the delay timer is active
            m = m - dt;                                     %then count it down by the time step
        end
    end

    
    if(mod(x,10) == 0)     %only store every tenth voltage to plot later
        V(n) = Vnew;
        S(n) = Snew;
        n = n + 1;
    end
end
 
%plot results
 figure(1); clf; plot(V, '-b'); hold on;

%initial conditions
V2(1) = Vresn;              %init cond for DSI cell 
n = 2;                     %keeps track of pointer in V and S vectors
m = tao;                   %default delay variable to tao seconds
delay = 1;                 %default initial delay trigger to on
Vnew = V2(1);               %start voltage at 0
Snew = 0;                  %start synaptic current at 0
omega = 3000;
epsilon = -.5;

%iterate through time steps to observe system behavior
for x=2:1:Last
    Vold = Vnew;           %keeps track of voltage from previous step  
    Sold = Snew;           %keeps track of old synaptic current from previous step
    
    if( x == omega)           %stimulate at time omega
       Vold = Vold + epsilon;
    end
    if ( delay == 1 && m <= 0)  %if the delay has been triggered and the tao-second countdown has expired
       Sold = 1;                    %then jump the synaptic current up to 1
       delay = 0;                   %and turn off the delay trigger
    end
    if(Vold >= Vthn)        %when the voltage hits threshold
        Vnew = 0;               %reset voltage to 0
        m = tao;                %start tao-second timer
        delay = 1;              %turn on delay trigger
    
    else                                            %if no other event happens, then simply
        Vnew = Vold + (-Vold + J + delV*Sold)*dt;       %change the voltage according to the DE
        Snew = Sold - (Sold/taoS)*dt;                   %change the synaptic current according to the exponential decay DE
        if(m>0)                                         %and if the delay timer is active
            m = m - dt;                                     %then count it down by the time step
        end
    end

    
    if(mod(x,10) == 0)     %only store every tenth voltage to plot later
        V2(n) = Vnew;
        S2(n) = Snew;
        n = n + 1;
    end
end
%plot results
 figure(1);  plot(V2, '-r');