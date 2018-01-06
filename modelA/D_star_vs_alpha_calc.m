%This code computes the costly performance of the remote-state estimation
%system with one transmitter and one remote estimator connected by an ideal
%channel
%Written by Jhelum Chakravorty
%Please see for theoretical reference the paper Jhelum Chakravorty and Aditya Mahajan, "Fundamental limits of remote estimation of autoregressive Markov processes 
%under communication constraints", IEEE Transactions on Automatic Control, March 2017.



clear all
clc

format long

p=0.3; %birth-rate

beta=[0.9;0.95;1.0]; %discount factor

color=['b';'r';'k'];

alpha=linspace(0.01,0.99,1000);

lambda=zeros(length(beta),1);
m=zeros(3,length(alpha));
j_star=zeros(3,1);
a=zeros(3,length(alpha));
D_opt=zeros(3,length(alpha));
theta=zeros(3,length(alpha));
D1=zeros(3,length(alpha));
D2=zeros(3,length(alpha));
N1=zeros(3,length(alpha));
N2=zeros(3,length(alpha));


for i=1:length(beta)
    
  for j=1:length(alpha)
        k=1;
        flag=0;
        L1=[];
        M1=[];
    while k>0 && flag==0
       if beta(i) ~= 1
          K=-2-(1-beta(i))/(beta(i)*p);
          lambda(i)=acosh(-K/2);
    
          [L1,M1]=calcLM(k,lambda(i),beta(i),p);
               
           a(i,j)=1/(1+alpha(j)-beta(i));
           if M1>= a(i,j)
               m(i,j)=k-1; %k_star in our writeup
               flag=1;
           end
        else M1=k^2/(2*p);
             a(i,j)=1/alpha(j);
             if M1 >= a(i,j)
                m(i,j)=k-1; %k_star in our writeup
                flag=1;
             end
       end
        k=k+1; 
    end
  end
end


%Calculation of D_star 
for i=1:length(beta)
for j=1:length(m)
    if m(i,j)~=0
     [L_1,M_1]=calcLM(m(i,j),lambda(i),beta(i),p);
     [L_2,M_2]=calcLM(m(i,j)+1,lambda(i),beta(i),p);
     N1(i,j)=1/M_1 - (1-beta(i)); N2(i,j) =1/M_2 -(1-beta(i));
     theta(i,j)=(alpha(j)-N2(i,j))/(N1(i,j)-N2(i,j));
     D1(i,j)=L_1/M_1;D2(i,j)=L_2/M_2;
     D_opt(i,j)=theta(i,j)*D1(i,j)+(1-theta(i,j))*D2(i,j);
    else D_opt(i,j)=0;
    end
end
end



for i=1:length(beta)
plot(alpha,D_opt(i,:),color(i));
hold on
end
ylim([0 3]);
set(gca,'PlotBoxAspectRatio',[5 3 1])
xlabel('$\alpha$','Interpreter','latex');
ylabel('$D^*_\beta(\alpha)$','Interpreter','latex');
legend('\beta = 0.9','\beta = 0.95','\beta = 0.1','Location','northeast','Interpreter','latex');     
