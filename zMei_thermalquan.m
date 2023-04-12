% define the comparable thermal quantity.
clear 
syms beta_m; % betal_m is modified beta which is defined 1/[(k*NA)*T]=1/(R*T), beta is 1/kT, k is the Boltzmann constant, and K 
% K is the Kelvin temperature. R = Na*k = 8.31. k =1.38*10^(-23). T =300;
% Na = 6.02*10^23. 
partionf = 0;
T = 300;
R = 8.31;
load('zv_sfc_2020.mat')
load('zz.mat')
[chang kuan]=size(V);
Thermalquan = 0;
for i=1:chang % determine the internal energy
    w = V(i,:)'*V(i,:);
    [evec,eval]=eig(w);
    max_eval(i) = max(diag(eval));
    %Thermalquan = Thermalquan + x(i)*max_eval;
end

for j =1:length(max_eval) % partition function z
    partionf = partionf + exp(-beta_m*max_eval(j));
end

 f(beta_m) = log(partionf);
 Energy = -diff(f(beta_m)); % Energy represents the internal energy;
 Entropy = R *(f - beta_m* diff(f(beta_m)));% Entropy
 FreeEnergy = - R*T*f(beta_m); % Helmholz free_engery
 
 
 beta_m =1/T*R; % beta = kT; k is the planck constant, and the T is the Kelvin temperature at the normal temperature.
  U = double(subs (Energy,beta_m))% subs is a function to calculate the Energy function by inputing beta; double
 % double is to show the results in fraction form. 
 S = double(subs(Entropy,beta_m))
 FE = double(subs(FreeEnergy,beta_m))
 
 se = 0;
 for k = 1: length(x)
     se =se +(-x(k)*log(x(k)));
 end