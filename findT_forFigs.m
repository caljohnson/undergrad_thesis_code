%findT_forFigs
function [T] = findT_forFigs(I, alpha, tau, tau_s)

fixedpt = 0;
%loop through values of delTk (spike time interval lengths)
delT = 0:.00001:10;
VF = zeros(size(delT,2),1);
VG = zeros(size(delT,2),1);
for i=1:size(delT,2)
    VF(i) = exp(delT(i)-tau) - I*(exp(delT(i)-tau)-1) + alpha*(tau_s/(tau_s-1))*(exp((delT(i)-tau)*(tau_s-1)/tau_s)-1);
    VG(i) = I*(1-exp(-tau)) - alpha*(tau_s)/(tau_s-1)*exp((-delT(i)+tau)/tau_s)*(exp(-tau/tau_s)-exp(-tau));
    
    %see if the two are close together
    if (VG(i)<VF(i)+.00025 && VG(i) > VF(i)-.00025)
        fixedpt = delT(i);
    end
    
end

    T = fixedpt;


end