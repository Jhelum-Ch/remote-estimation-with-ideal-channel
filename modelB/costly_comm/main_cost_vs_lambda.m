%This code computes the costly performance of the remote-state estimation
%system with one transmitter and one remote estimator connected by an ideal
%channel
%Written by Jhelum Chakravorty
%Please see for theoretical reference the paper Jhelum Chakravorty and Aditya Mahajan, "Fundamental limits of remote estimation of autoregressive Markov processes 
%under communication constraints", IEEE Transactions on Automatic Control, March 2017.



clear all
clc
format long

var  = 1;
betaVec = [0.9;0.95;1];
a=1; %the parameter 'a' in the dynamics

lambda_vec=[50,100,150,200,250,300,350,400,450,500,550,600,650,700];
C=zeros(length(betaVec),length(lambda_vec));
tic

   
for i=1:length(betaVec)
      beta = betaVec(i);
  [C(i,:),k_opt_vec] = compute_cost(lambda_vec, beta, var,a);
  

set(gca,'PlotBoxAspectRatio',[5 3 1])
% For plotting through ssh
% set(fig, 'Visible', 'off')

color=['b';'r';'k'];

plot(lambda_vec, C(i,:),color(i));
hold on
end
toc

xlabel('$\lambda$','Interpreter','latex');
ylabel(' $C^*_\beta(\lambda)$','Interpreter','latex');
legend('\beta = 0.9','\beta = 0.95','\beta = 0.1','Location','northwest','Interpreter','latex');     
