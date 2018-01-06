%This code computes the costly performance of the remote-state estimation
%system with one transmitter and one remote estimator connected by an ideal
%channel
%Written by Jhelum Chakravorty
%Please see for theoretical reference the paper Jhelum Chakravorty and Aditya Mahajan, "Fundamental limits of remote estimation of autoregressive Markov processes 
%under communication constraints", IEEE Transactions on Automatic Control, March 2017.



clear all
clc

format long

p=0.3; % birth-rate

beta=[0.9;0.95;1.0]; %discount factor

c_final=25; %SEt an upper limit of transmission cost

lambda=zeros(length(beta),1);
 
for i=1:length(beta)
    D=-2-(1-beta(i))/(beta(i)*p);
         lambda(i)=acosh(-D/2);

        k=1; %thresholds
        c=0.1;
    if beta(i)==1
        lambda_beta=[]; %initialize vector of critical Lagrange parameter; 
        L=[];
        M=[];
     while k>0 && c<c_final  
     lambda_beta(k)=(k*(k+1)*(k^2+k+1))/(6*p*(2*k+1));
     [L(i,k),M(i,k)]=calcLM(k,lambda(i),beta(i),p);
     c=lambda_beta(end);
     C(k) = (L(i,k)+lambda_beta(k))/M(i,k);
     k=k+1;
     end
     lambda_beta1=lambda_beta;
    elseif beta(i)==0.9
        lambda_beta=[];
        L=[];
        M=[];
        

         
         while k>0 && c<c_final

         [L(i,k),M(i,k)]=calcLM(k,lambda(i),beta(i),p);
         D(i,k)=L(i,k)/M(i,k);
         if k>1 && c<c_final
         lambda_beta(k-1)=(D(i,k) - D(i,k-1))/(1/M(i,k-1) - 1/M(i,k));
         c=lambda_beta(k-1);
         C_disc09(k-1) = L(i,k-1)/M(i,k-1)+ c*((1/M(i,k-1))-(1-beta(i)));
         lambda_beta09(k-1)=lambda_beta(k-1); 
         end
            k=k+1;
         end
    else  lambda_beta=[];
        L=[];
        M=[];
        

         
         while k>0 && c<c_final
             [L(i,k),M(i,k)]=calcLM(k,lambda(i),beta(i),p);
             D(i,k)=L(i,k)/M(i,k);
         if k>1 && c<c_final
         lambda_beta(k-1)=(D(i,k) - D(i,k-1))/(1/M(i,k-1) - 1/M(i,k));
         c=lambda_beta(k-1);
         C_disc095(k-1) = L(i,k-1)/M(i,k-1)+ c*((1/M(i,k-1))-(1-beta(i)));
         lambda_beta095(k-1)=lambda_beta(k-1); 
         end       
            k=k+1;
         end
    end

end




%Plots
lambda_tick=[1, 5, 10, 15, 20, 25, 30, 35];

plot(lambda_beta09,C_disc09, 'b')
hold on
plot(lambda_beta095,C_disc095,'r')
hold on
plot(lambda_beta1,C,'k')

set(gca, 'XTick', lambda_tick); % Change x-axis ticks
set(gca, 'XTickLabel', lambda_tick); 
set(gca,'PlotBoxAspectRatio',[5 3 1])
xlabel('$\lambda$','Interpreter','latex');
ylabel(' $C^*_\beta(\lambda)$','Interpreter','latex');
legend('\beta = 0.9','\beta = 0.95','\beta = 0.1','Location','southeast','Interpreter','latex');