%file neural_coupling_thesis_part1.m
%Carter Johnson, Undergrad Thesis
%Coupling in neural populations using the
%Morris-Lecar model of excitable barnacle muscle fiber
%adapted from Morris and Lecar (1981) Biophysical Journal 35 pp. 193-231
%Figures 8.6 and 8.7

function neural_coupling_thesis_part1
% PARAMETERS  
%declare model parameters
global C;
global gbarCa;
global ECa;
global gbarK;
global EK;
global gleak;
global Eleak;
global v1;
global v2;
global v3;
global v4;
global phi;
global tau;
global Iapplied;
%parameter values
C= 20 ; %microfarad/cm^2 
gbarCa= 4.4; % millisiemens/ cm^2 
ECa= 120; %millivolts
gbarK= 8;% millisiemens/ cm^2 
EK= -84; %millivolts
gleak= 2;% millisiemens/ cm^2 
Eleak= -60;%millivolts
v1= -1.2; %millivolts
v2= 18 ; %millivolts
v3= 2 ; %millivolts
v4= 30; %millivolts
phi = 0.04; % per millisecond
tau=0.8;
Iapplied=150;

%% Figure 1 
%FIND SOLUTION TRAJECTORY

%generate nullclines
figure(1)
hold on
myfun1 = @(Vn,wn) (1/C)*(gbarCa*(0.5*(1+tanh((Vn-v1)/v2)))*(ECa-Vn) + gbarK*wn*(EK-Vn) + gleak*(Eleak-Vn)+Iapplied);
%myfun2 = @(Vn,wn) phi*((0.5*(1+tanh((Vn-v3)/v4)))-wn)/((1/(cosh((Vn-v3)/(2*v4)))));
myfun2 = @(Vn,wn) phi*((0.5*(1+tanh((Vn-v3)/v4)))-wn)/tau;
a1=ezplot(@(Vn,wn) myfun1(Vn,wn), [-80 60 -0.2 1]);
setcurve('color','black','linestyle','--', 'Linewidth', 1) 
set(gca, 'fontsize', 14)
a2=ezplot(@(Vn,wn) myfun2(Vn,wn), [-80 60 -0.2 1]);
setcurve('color',[0.5 0.5 0.5], 'Linewidth', 1) 
legend('V nullcline', 'w nullcline')
title('')

%USE Euler's method to calculate trajectory
%Euler's method constants
dt = .1;         %time step
t_last = 10000;     %max steps
X = zeros(t_last,2);  %initialize trajectory vector
X(1,:) = [-40, 0];  %initial conditions
for i = 2:t_last
    dX = morris_lecar_ddt(i,X(i-1,:));  %get derivative
    X(i,1) = dX(1)*dt + X(i-1,1);       %next voltage V
    X(i,2) = dX(2)*dt + X(i-1,2);       %next recovery var W
end
%plot trajectory in phase plane on top of nullclines
figure(1); plot(X(:,1), X(:,2), 'k', 'Linewidth', 2); 
xlabel('Membrane Voltage V', 'fontsize',12)
ylabel('Potassium gating variable W', 'fontsize',12), axis([-70 60 0. 0.75])


%% Figures 2 and 3 
%FIND LIMIT CYCLE
nmax = 1000; %max number of time steps in one period
vmark = -20.0; %arbitrary init. and final V, used to check LC
x0(1,:) = X(1,:);  %initial state X=(V,w)
dx=10000.0;   %arbitrarily large starting "precision", need to start loop
dxcrit=.1;    %precision for periodic orbit
i_LC = 1;     %starting point index of LC search
while (dx>dxcrit) 
% iterate until periodic orbit found to desired precision (dxcrit).
    x=x0;                   %initialize starting point of LC
    ii=1;                   %start counter of total # of points in LC
    flagg=0;                %have not found LC endpoint yet   
    while (flagg==0 && ii<=nmax) % time-step until system returns to
                                 % v=vmark with dv/dt>0                               
        y=X(i_LC+ii,:);             %get next step in trajectory     
        if (x(1)<vmark && y(1)>=vmark)  %check if we return to Vmark
            flagg=1;   %exit while loop since we may have found LC endpoint
        end;
        x=y;                             %copy LC point
        X_LC(ii,:)=x;                    %output LC point
        ii=ii+1;                         %up total # of points in LC
    end; 
    dx = norm(x-x0);                 %check how close our endpoints are
    x0=x;                            %copy endpoint
    i_LC = ii+i_LC;                  %set index of LC search
end;

% rearrange output so that minimum v is at t=0.
ii=ii-1;
[~,imin]=min(X_LC(:,1));
X_LC_2(1:ii-imin+1,:)=X_LC(imin:ii,:);
X_LC_2(ii-imin+2:ii,:)=X_LC(1:imin-1,:);
X_LC=X_LC_2;
T= ii;      %get vector length of LC

% OUTPUT:   X_LC(1:ii,:) contains periodic orbit
%           (ii-1)*dt=~period of orbit
figure(2); clf; plot(X_LC(:,1), 'r','Linewidth', 2); 
ylabel('V', 'fontsize',12), xlabel('t', 'fontsize',12), axis([0 T -50 40])
figure(3); clf; plot(X_LC(:,2), 'b','Linewidth', 2); 
ylabel('W', 'fontsize',12), xlabel('t', 'fontsize',12), axis([0 T 0 .65])

%% Figure 4 
%FIND THE PRC
%use limit cycle to solve for the PRC by
%solving the adjoint equation
%dZ/dt = -J^t|(X_LC) * Z

%set up PRC variable
Z = zeros(T,2);
Z(T,:) = [-.2, 50];
    %loop through solution backwards multiple times to ensure convergence
    %to limit cycle
for j = 1:100
    %solve adjoint equation backwards in time with Euler's method
for i = T:-1:1
    v = X_LC(i,1);
    w = X_LC(i,2);  
    %create the transposed Jacobian matrix at the specific point on the
    %limit cycle of (V,W)
    Fv = 1/C*( -gbarCa * 0.5*(1+tanh((v - v1)/v2)) + gbarCa*(1/(2*v2))*(sech((v - v1)/v2))^2*(ECa - v) - gbarK * w - gleak);
    Fw = 1/C * gbarK * (EK - v);
    Gv = -phi/tau * (sech((v-v3)/v4))^2*(1/(2*v4));
    Gw = phi/tau;

    %calculate next value of Z
    if i > 1   
     Z(i-1,1) = -(Fv*Z(i,1) - Gv*Z(i,2))*(-dt) + Z(i,1);
     Z(i-1,2) = -(Fw*Z(i,1) - Gw*Z(i,2))*(-dt) + Z(i,2);
    else
     Z(T,1) = -(Fv*Z(1,1) - Gv*Z(1,2))*(-dt) + Z(1,1);
     Z(T,2) = -(Fw*Z(1,1) - Gw*Z(1,2))*(-dt) + Z(1,2);   
    end
    
end
end
figure(4); plot(Z(:,1));

%normalize PRC
%calculate X'_LC, derivative of limit cycle
dX_LC = diff(X_LC(1:T,:))/dt;
%match derivatives at endpoints
dX_LC = [dX_LC; dX_LC(1,:)];
%take average of the inner product of X'_LC and Z
z0 = zeros(T,1);
for i = 1:T
z0(i) =  dX_LC(i,:)*transpose(Z(i,:));
end
aveZ0 = median(z0);
%normalize PRC by aveZ0
Z = Z/aveZ0;

%plot PRC
figure(4); subplot(2,1,1); plot(Z(:,1), 'r', 'Linewidth', 2);
ylabel('V', 'fontsize',12), xlabel('t', 'fontsize',12), axis([0 T -.5 1])
subplot(2,1,2); plot(Z(:,2), 'b', 'Linewidth', 2);
ylabel('W', 'fontsize',12), xlabel('t', 'fontsize',12), axis([0 T -100 100])

%% Figure 5
%CALCULATE phase-locked states of neuron with forcing field

%set forcing electrical field Ie(t)
Amp = -600;  %amplitude of forcing field
Ie = zeros(T,1);
for i = 1:T
Ie(i) = Amp*sin(2*pi*i/T);
end

 %copy Ie to do convolution
 Ie_copy = Ie;
 %set phase model H(t);
 H = Z(:,1);
 
%convolute iPRC with interaction function
for psi = 1:T
  %shift Ie by -psi+1
    %start Ie at phase psi-1
    for i = psi:T
     Ie(i) = Ie_copy(i-psi+1);
    end
    %end S at phase psi-1
    for i = 1:psi-1
     Ie(i) = Ie_copy(T-(psi-1)+i);
    end
    
    %compute convolution integral
    intSum = Z(1,1)*Fv*Ie(1);
    for tHat = 2:T
        intSum = intSum + Z(tHat,1)*Fv*Ie(tHat);
    end
    H(psi) = intSum/T;
end
 
 %calculate dphi/dy = H(-phi)-H(phi))
 %dphi= -H+H2;
 figure(5); clf; plot(H, 'k', 'LineWidth', 2);
 ylabel('d \Phi / dt', 'fontsize',30), xlabel('t', 'fontsize',30)
 
 %% Figure 6
 %FIND phase-locked position and then
 %SHOW phase-locked limit cycle timings
 %search dPhi/dt for zero
 old_D = H(1);
 for i = 2:T
     new_D = H(i);
     %check if zero is crossed
     if (old_D > 0 && new_D < 0 || new_D == 0)
         %choose the closer index
         if (abs(old_D) > abs(new_D))
           SS = i;  %set index of zero/stable state
           break;
         else
           SS = i-1;
           break;
         end
     end
     old_D = new_D;
 end
 
 %reindex limit cycle to time with input current
 retimedV = zeros(T,1); 
 for i = 1:T-SS
     retimedV(i) = X_LC(i+SS,1);
 end
 for i = T-SS+1 : T
     retimedV(i) = X_LC(i - T + SS , 1);
 end
 
 %show timings
 figure(6); clf; subplot(2,1,1); plot( retimedV, 'r', 'LineWidth', 2);
 ylabel('V', 'fontsize',30), axis([0 T -50 50]),
 subplot(2,1,2); plot(Ie_copy, 'b', 'LineWidth', 2);
 ylabel('I', 'fontsize',30), xlabel('t', 'fontsize',30), axis([0 T -50 50]),
 
 %% Figure 7
 %CALCULATE phase-locking in morris-lecar pairs 
 %with exponentially-decaying synaptic coupling
  
%set function S(t)
 S = zeros(T,1);
 for i = 1:T
     S(i) = exp(-(i-1)/100);
 end
 S_copy = S;
 
 %set function H(t)
 H3 = zeros(T,1);
%convolute iPRC with synaptic input
for psi = 1:T
  %shift S by -psi+1
    %start S at phase psi-1
    for i = psi:T
     S(i) = S_copy(i-psi+1);
    end
    %end S at phase psi
    for i = 1:psi-1
     S(i) = S_copy(T-(psi-1)+i);
    end
    
    %compute convolution integral
    intSum = Z(1)*eps*S(1)/2 +  Z(T)*eps*S(T)/2;
    for tHat = 2:T-1
        intSum = intSum + Z(tHat)*eps*S(tHat);
    end
    H3(psi) = intSum/T;
end
 
 %calc H2 = H3(-t)
 H2 = zeros(T,1);
 for i = 1:T
     H2(i) = H3(T+1-i,1);
 end

 dphi2= H3-H2;
 figure(7); clf; plot(dphi2, 'k', 'LineWidth', 2);
 ylabel('d \Phi / dt', 'fontsize',30), xlabel('t', 'fontsize',30);
 %% Figure 8
 %FIND phase-locked position and then
 %SHOW phase-locked limit cycle timings
 
 %search dPhi/dt for zero
 old_D = dphi2(1);
 for i = 2:T
     new_D = dphi2(i);
     %check if zero is crossed
     if (old_D > 0 && new_D < 0 || new_D == 0)
         %choose the closer index
         if (abs(old_D) > abs(new_D))
           SS = i;  %set index of zero/stable state
           break;
         else
           SS = i-1;
           break;
         end
     end
     old_D = new_D;
 end
 
 %reindex limit cycle to time with input current
 retimedV = zeros(T,1); 
 for i = 1:T-SS
     retimedV(i) = X_LC(i+SS,1);
 end
 for i = T-SS+1 : T
     retimedV(i) = X_LC(i - T + SS , 1);
 end
 
 %show timings
 figure(8); clf; subplot(2,1,1); plot( retimedV, 'r', 'LineWidth', 2);
 ylabel('V1', 'fontsize',30), axis([0 T -50 50]),
 subplot(2,1,2); plot(X_LC(:,1), 'b', 'LineWidth', 2);
 ylabel('V2', 'fontsize',30), axis([0 T -50 50]), xlabel('t', 'fontsize',30)

  %% Figures 9-13
 %DEMONSTRATE phase-locking in morris-lecar groups of 4 
 %with exponentially-decaying synaptic coupling
  
 %show timings
 figure(10); clf; subplot(4,1,2); plot( retimedV, 'r', 'LineWidth', 2);
 ylabel('V2', 'fontsize',30), axis([0 T -50 50]),
 subplot(4,1,1); plot(X_LC(:,1), 'b', 'LineWidth', 2);
 ylabel('V1', 'fontsize',30), axis([0 T -50 50]),
 subplot(4,1,3);plot( retimedV, 'r', 'LineWidth', 2);
 ylabel('V3', 'fontsize',30), axis([0 T -50 50]),
 subplot(4,1,4); plot(X_LC(:,1), 'b', 'LineWidth', 2);
 ylabel('V4', 'fontsize',30), axis([0 T -50 50]),xlabel('t', 'fontsize',30)
 
 figure(11); clf; subplot(4,1,4); plot( retimedV, 'r', 'LineWidth', 2);
 ylabel('V4', 'fontsize',30), axis([0 T -50 50]), xlabel('t', 'fontsize',30)
 subplot(4,1,2); plot(X_LC(:,1), 'b', 'LineWidth', 2);
 ylabel('V2', 'fontsize',30), axis([0 T -50 50]),
 subplot(4,1,3);plot( retimedV, 'r', 'LineWidth', 2);
 ylabel('V3', 'fontsize',30), axis([0 T -50 50]),
 subplot(4,1,1); plot(X_LC(:,1), 'b', 'LineWidth', 2);
 ylabel('V1', 'fontsize',30), axis([0 T -50 50]),
 
 figure(12); clf; subplot(4,1,1); plot( X_LC(:,1), 'b', 'LineWidth', 2);
 ylabel('V1', 'fontsize',30), axis([0 T -50 50]),
 subplot(4,1,2); plot(X_LC(:,1), 'b', 'LineWidth', 2);
 ylabel('V2', 'fontsize',30), axis([0 T -50 50]),
 subplot(4,1,3);plot( X_LC(:,1), 'b', 'LineWidth', 2);
 ylabel('V3', 'fontsize',30), axis([0 T -50 50]),
 subplot(4,1,4); plot(X_LC(:,1), 'b', 'LineWidth', 2);
 ylabel('V4', 'fontsize',30), axis([0 T -50 50]), xlabel('t', 'fontsize',30)
 
 %%%%%% 4 cell 1/4 phase synchrony
 phi_diff = T/4;
 %reindex limit cycles to time 
 xlc2 = zeros(T,1); 
 for i = 1:T-phi_diff
     xlc2(i) = X_LC(i+phi_diff,1);
 end
 for i = T-phi_diff+1 : T
     xlc2(i) = X_LC(i - T + phi_diff , 1);
 end
 
  xlc3 = zeros(T,1); 
 for i = 1:T-2*phi_diff
     xlc3(i) = X_LC(i+2*phi_diff,1);
 end
 for i = T-2*phi_diff+1 : T
     xlc3(i) = X_LC(i - T + 2*phi_diff , 1);
 end
 xlc4 = zeros(T,1); 
 for i = 1:T-3*phi_diff
     xlc4(i) = X_LC(i+3*phi_diff,1);
 end
 for i = T-3*phi_diff+1 : T
     xlc4(i) = X_LC(i - T + 3*phi_diff , 1);
 end
 
  figure(13); clf; subplot(4,1,1); plot( X_LC(:,1), 'r', 'LineWidth', 2);
 ylabel('V1', 'fontsize',30), axis([0 T -50 50]),
 subplot(4,1,2); plot(xlc4(:,1), 'r', 'LineWidth', 2);
 ylabel('V2', 'fontsize',30), axis([0 T -50 50]),
 subplot(4,1,3);plot( xlc3(:,1), 'r', 'LineWidth', 2);
 ylabel('V3', 'fontsize',30), axis([0 T -50 50]),
 subplot(4,1,4); plot(xlc2(:,1), 'r', 'LineWidth', 2);
 ylabel('V4', 'fontsize',30), axis([0 T -50 50]), xlabel('t', 'fontsize',30)
end

%dynamics
function dS = morris_lecar_ddt(t,S)

global C;
global gbarCa;
global ECa;
global gbarK;
global EK;
global gleak;
global Eleak;
global v1;
global v2;
global v3;
global v4;
global phi;
global tau;
global Iapplied;

%locally define state variables:
V=S(1);
w=S(2);

%local functions:
m_inf = (0.5*(1+tanh((V-v1)/v2)));
w_inf = (0.5*(1+tanh((V-v3)/v4)));

ddt_V = (1/C)*(gbarCa*m_inf*(ECa-V) + gbarK*w*(EK-V) + gleak*(Eleak-V)+Iapplied);
ddt_w = phi*(w_inf-w)/(tau);

dS=[ddt_V, ddt_w];

end

%change properties of last curve in current figure
function setcurve(varargin)
h=get(gca,'children');
set(h(1),varargin{:})
end