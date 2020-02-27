clc
close all

snr = 40;
lambda0 = 100;
T = 30*9*(10^-6);
d = 500*(10^-9);
r = 45*(10^-9);
L=5;
time_slot = [1:1000];
capD = 4.265*(10^-10); % diffusion coefficient
i=1;
ri = zeros(1,1000);
sj = randn(1, length(time_slot))>0.5;

P_0 = (r/d)*(erfc((d-r)/sqrt(4*capD*i*T))-erfc((d-r)/sqrt(4*capD*(i-1)*T))); %constant
Ntx = 2.*lambda0.*T.*(10.^(snr./10))./P_0;


for i = 1:length(time_slot)
	P_i1(i) = (r/d)*(erfc((d-r)/sqrt(4*capD*i*T))-erfc((d-r)/sqrt(4*capD*(i-1)*T)));
end

for j = 1:length(time_slot)
	cj(j) = Ntx.*P_i1(j);
end

c0 = 54;
sum1 = sum(sj.*cj);
avg = (lambda0*T + sum1);
avg1 = avg + sj.*c0;
ri(1,:) = poissrnd(avg1);

   %tau(k) = cj(1)/log(1 + (cj(1)/(sum(cj/2)+(lambda0*T))));


%ri = reshape(ri, [1, 56*5])

%ri = reshape
% sj = repmat(sj, 1, 56)

A = [ri; sj];
writematrix(A, 'MC_Ntx_40.csv');

Ntx
ri
sj


%writematrix(sj, 'MC_Ntx.csv');