%This code computes the costly performance of the remote-state estimation
%system with one transmitter and one remote estimator connected by an ideal
%channel
%Written by Jhelum Chakravorty
%Please see for theoretical reference the paper Jhelum Chakravorty and Aditya Mahajan, "Fundamental limits of remote estimation of autoregressive Markov processes 
%under communication constraints", IEEE Transactions on Automatic Control, March 2017.




clear all
clc
format long

beta=[0.9;0.95;1];  %values of discount factor
tol=10^(-2);

T_vec=1:1:10;
alpha1_vec=1./T_vec;

alpha2_vec=zeros(1,8);
for i=1:8
    alpha2_vec(i)=(i+1)/(i+2);
end


alpha_rev=fliplr(alpha1_vec(1,2:end));

alpha2_vec=[alpha2_vec 1];

alpha = [linspace(0.01,0.1,9) linspace(0.15, 0.9, 16)];

var=1;



RHS=@(e) e.^0;
RHS1=@(e) e.^2;



AbsTol = 1e-6;
RelTol = 1e-3;

color=['b';'r';'k'];



k_star=zeros(length(var),length(alpha));
D_k_star=zeros(length(var),length(alpha));


for j=1:length(beta)
    M=[];
t=[];
L=[];
    
    
kernel=@(e,x) (1/sqrt(2*pi*var))*exp(-(x-e).^2/(2*var)); % Gaussian

%For uniform distribution
state_max=15;

for i=1:length(alpha)
    

    k_min=0; k_max=100;
    
    keepLooping = true;
    

while keepLooping 

k_guess = 0.5*(k_min+k_max);

%Solve M for k_guess
a=-k_guess; b=k_guess;
scalar=1/beta(j);
rescaled_RHS=@(e) RHS(e)/beta(j);
[sol,errest,cond] = Fie(scalar,a,b,1,kernel,rescaled_RHS,AbsTol,RelTol);
M=sol.x;
t = sol.s;
M_k_guess=M((length(M)+1)/2);

if abs(M_k_guess - 1/(alpha(i)-beta(j)+1)) < tol
  k_star(i) = k_guess;
  break; %to break out of the while loop
elseif M_k_guess < 1/(alpha(i)-beta(j)+1)
   k_min = k_guess;
 else
   k_max = k_guess;
end
end
a=-k_star(i); b=k_star(i);
scalar=1/beta(j);
rescaled_RHS1=@(e) RHS1(e)/beta(j);
[sol,errest,cond] = Fie(scalar,a,b,1,kernel,rescaled_RHS1,AbsTol,RelTol);
L=sol.x;
t = sol.s;
L_k_star(j,i)=L((length(L)+1)/2);

D_k_star(j,i)=L_k_star(j,i)*(alpha(i)-beta(j)+1);  %D-T function
end
end


 for j=1:length(beta)
     plot(alpha,D_k_star(j,:),color(j));
      hold on
 end
 
     set(gca,'PlotBoxAspectRatio',[5 3 1])
     xlabel('$\alpha$','Interpreter','latex');
     ylabel('$D^*_\beta(\alpha)$','Interpreter','latex');
     legend('\beta = 0.9','\beta = 0.95','\beta = 0.1','Location','northeast','Interpreter','latex');     



