%%Figures 2 and 3
    %generates figures 2 and 3, the evolution of the standard LIF and DSI cells
    
%parameters

%bias current I
    I = 1.1;
%stable fixed period length T
    T = 2;
%synaptic strength alpha (-delV)
   delV = .5;
%synaptic delay tao
    tao = 1.1;
%decay rate of synaptic current taoS
    taoS = 2;
   
%Evolution of Two Models side-by-side    
    
%Standard LIF Model stuff
T2 = -log(1-1/I); %period

%simulation constants
dt = 0.001; %time step
tF = 10.5*T2;%final time
Last = floor(tF/dt + 1);

%standard LIF
V1 = zeros(Last,1); %create space for V
V2 = V1;
S2 = V1;
t = 0:dt:dt*floor(Last - 1);

%DSI cell
n = 2;                          %keeps track of pointer in V and S vectors
m = tao;                        %default delay variable to tao seconds
delay = 1;                      %default initial delay trigger to on
Vold = V2(1);                    %start voltage at 0
Sold = S2(1);                    %start S at steady state, value at Tstar

for x=2:1:Last
    
    if(V1(x-1)>=1)
     V1(x) = 0;
    else
     V1(x)= V1(x-1) + (-V1(x-1) + I)*dt;
    end
    
    Vnew = Vold + (-Vold + I - delV*Sold)*dt;       %change the voltage according to the DE
    Snew = Sold - (Sold/taoS)*dt;                   %change the synaptic current according to the DE
    if(m>0)                                         %and if the delay timer is active
      m = m - dt;                                     %then count it down by the time step
    end   
    
    if ( delay==1 && m<= 0)  %if the delay has been triggered and the tao-second countdown has expired
       Snew = 1;                    %then jump the synaptic current up to 1
       delay = 0;                       %and turn off the delay trigger
    end
  
    
    if(Vnew >= 1)                          %when the voltage hits threshold
        m = tao;                              %start tao-second timer
        delay = 1;                            %turn on delay trigger
        dtstar = dt*(1-Vold)/(Vnew-Vold);  
        Vnew = 0;                             %reset voltage to 0 
        Snew = Sold - (Sold/taoS)*dtstar;     %change the synaptic current according to the DE
           
    end
    
    Vold = Vnew;                %keeps track of voltage from previous step  
    Sold = Snew;                %keeps track of old synaptic current from previous step
    
    V2(x) = Vnew;
    S2(x) = Snew;
        
end


%PLOT DETAILS

 figure(2); subplot(3,1,1); LIF = plot(t,V1, '-k'); set(LIF, 'LineWidth', 2); title('(a)'); axis([0,25,0,1.1]); xlabel('time'); ylabel('V');
 subplot(3,1,2); DSI = plot(t,V2, '-b'); set(DSI, 'LineWidth', 2); title('(b)'); axis([0,25,0,1.1]); xlabel('time'); ylabel('V');
 subplot(3,1,3); S = plot(t, S2, '-r'); set(S, 'LineWidth', 2); axis([0,25,0,1.1]); title('(c)'); xlabel('time'); ylabel('S');
 
