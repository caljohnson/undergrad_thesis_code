%Period dependence
    %the purpose here is to show how the stable period depends on the
    %various parameters- applied current I, delay timing tao, synapse
    %weight alpha, synaptic time constant tao_s
    
%simulation constants    
I = 1.1;                %applied current
Tnat = log(I/(I-1));    %natural period given I
dt = 0.1;              %time step of delay timing

%loop through delay timings

    %delay timing loop's constants
    taoS = 2;
    tao = zeros(Tnat/dt, 1);
    alpha = zeros(1/dt + 1, 1);
   
    T = zeros(Tnat/dt, 1/dt+1); 
    SC = zeros(Tnat/dt, 1/dt+1);  %stabilityCheck G'/F'(T*)
  
for i = 1:Tnat/dt
    tao(i) = i*dt;
    Vtao = I*(1-exp(-tao(i)));
    for j = 1:1/dt+1
        alpha(j) = (j-1)*dt;
        [T(i, j),SC(i,j)] = periodLogMap(I, -alpha(j), tao(i), taoS, 0);
    end
end

figure(1); clf; surf(alpha, tao, T);
figure(2); clf; surf(alpha, tao, abs(SC));
