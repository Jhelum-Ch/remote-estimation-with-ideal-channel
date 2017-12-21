clear all 
clc


Num_of_samples = 10000;

flag_dis = true;

p = 0.3; 
 
beta = 1.0; 
epsilon_vec=0;
%epsilon_vec = linspace(0,0.9,10);
color = ['b','m','c','g','r','k'];

K=8; %range of threshold values
K_vec = linspace(1,K,K);

%lambda = [20, 20]


for i=1:length(epsilon_vec)
    epsilon=epsilon_vec(i);
    %lambda = zeros(1,length(K_vec));

for k=1:length(K_vec)
M_sum = 0;
L_sum = 0;
N_sum = 0;
ITER = 10000;
for count = 1:Num_of_samples
  if flag_dis
     [n,t, D, iter] = sampleM_dis_MC(beta,p,k,epsilon,ITER);
  else
     [n,t, D] = sampleM_MC(p,k);
  end
  M_sum = M_sum + t;
  L_sum = L_sum + D;
  N_sum = N_sum + n;
end


M_est(i,k) = M_sum/ Num_of_samples;
L_est(i,k) = L_sum/ Num_of_samples;
%N_est(i,k) = (1-beta)*N_sum/ (Num_of_samples);
%N_test(i,k) = ((1-beta^(1/(1-epsilon)))/(1-beta))*(1/M_est(i,k)-(1-beta));
N_test(i,k) = 1/M_est(i,k);
end
end
