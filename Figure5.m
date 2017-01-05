%Figure5

%NEWFEQGMAP We plot F and G on the same axis and look for points of
%intersection (the fixed point)

I = 1.1;
alpha = .5;
tau = 1.5;
tau_s = 2;

fixedpt = 0;

%loop through values of delTk (spike time interval lengths)
delT = 0:.001:10;
VF = zeros(size(delT,2),1);
VG = zeros(size(delT,2),1);
for i=1:size(delT,2)
    VF(i) = exp(delT(i)-tau) - I*(exp(delT(i)-tau)-1) + alpha*(tau_s/(tau_s-1))*(exp((delT(i)-tau)*(tau_s-1)/tau_s)-1);
    VG(i) = I*(1-exp(-tau)) - alpha*(tau_s)/(tau_s-1)*exp((-delT(i)+tau)/tau_s)*(exp(-tau/tau_s)-exp(-tau));
    if VF(i)==VG(i)
        fixedpt = delT(i);
    end
end

figure(5); plot(delT, VF, 'r--'); hold on; plot(delT, VG, 'b');
axis([2,7,0,3]);
xlabel('\Delta T'); ylabel('V'); title('F(\Delta T) vs G(\Delta T)');
legend('F(\Delta T)', 'G(\Delta T)', 'Location', 'northwest');

