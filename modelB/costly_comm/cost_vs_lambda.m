clear all
clc
format long

var  = 1;
%beta = [0.9;1];
beta = 0.9;
%A=0.8;
%A=[0.8;1;2]; %values of 'a'

A=1;
lambda_vec=50;

%lambda_vec=[50,100,150,200,250,300,350,400,450,500,550,600,650,700];
%lambda_vec= [0.1, 0.5, 1, 5, 10, 20, 50, 100];
C=zeros(length(A),length(lambda_vec));
tic
for i=1:length(A)
    a=A(i);
   
   %for i=1:length(beta)
[C(i,:),k_opt_vec] = compute_cost(lambda_vec, beta, var,a);
   %end
end
toc


% fig = figure();
% 
% set(gca,'PlotBoxAspectRatio',[5 3 1])
% % For plotting through ssh
% % set(fig, 'Visible', 'off')
% 
% color=['b';'r';'k'];
% 
% 
% for i=1:length(A)
% plot(lambda_vec, C(i,:),color(i));
% % for i=1:length(beta)
%  hold on
% % end
% end
% xlabel('$\lambda$','Interpreter','latex');
% ylabel(' $C^*_\beta(\lambda)$','Interpreter','latex');
%  
% print(fig,'PlotofC_opt','-dpdf')    
