clear all
clc
format long

lambda=0.9; %beta
%lambda=[0.9;0.95;1];  %value of beta
tol=10^(-2);

T_vec=1:1:10;
alpha1_vec=1./T_vec;

alpha2_vec=zeros(1,8);
for i=1:8
    alpha2_vec(i)=(i+1)/(i+2);
end


alpha_rev=fliplr(alpha1_vec(1,2:end));

alpha2_vec=[alpha2_vec 1];

%alpha=[alpha_rev alpha2_vec];
%alpha = linspace(0.01, 0.9,40);
alpha = [linspace(0.01,0.1,9) linspace(0.15, 0.9, 16)];

%alpha =  linspace(0.15, 0.95, 17);

%p=0.3;

%var=2*p;
%var=[0.5;1.0;10];
var=1;



RHS=@(e) e.^0;
RHS1=@(e) e.^2;



AbsTol = 1e-6;
RelTol = 1e-3;

color=['b';'r';'k'];



k_star=zeros(length(var),length(alpha));
D_k_star=zeros(length(var),length(alpha));


for j=1:length(lambda)
    M=[];
t=[];
L=[];
    
    
kernel=@(e,x) (1/sqrt(2*pi*var))*exp(-(x-e).^2/(2*var)); % Gaussian

%For uniform distribution
state_max=15;

%kernel=@(e,x) (abs((x-e)) < state_max)*(1/(2*state_max)); % Uniform

for i=1:length(alpha)
    

    k_min=0; k_max=100;
    
    keepLooping = true;
    

while keepLooping 

k_guess = 0.5*(k_min+k_max);

%Solve M for k_guess
a=-k_guess; b=k_guess;
%scalar=lambda(j);
scalar=1/lambda(j);
rescaled_RHS=@(e) RHS(e)/lambda(j);
%[sol,errest,cond] = Fie(scalar,a,b,1,kernel,RHS,AbsTol,RelTol);
[sol,errest,cond] = Fie(scalar,a,b,1,kernel,rescaled_RHS,AbsTol,RelTol);
M=sol.x;
t = sol.s;
M_k_guess=M((length(M)+1)/2);

if abs(M_k_guess - 1/(alpha(i)-lambda(j)+1)) < tol
  k_star(i) = k_guess;
  break; %to break out of the while loop
elseif M_k_guess < 1/(alpha(i)-lambda(j)+1)
   k_min = k_guess;
 else
   k_max = k_guess;
end
end
a=-k_star(i); b=k_star(i);
%scalar=lambda(j);
scalar=1/lambda(j);
rescaled_RHS1=@(e) RHS1(e)/lambda(j);
%[sol,errest,cond] = Fie(scalar,a,b,1,kernel,RHS1,AbsTol,RelTol);
[sol,errest,cond] = Fie(scalar,a,b,1,kernel,rescaled_RHS1,AbsTol,RelTol);
L=sol.x;
t = sol.s;
L_k_star(j,i)=L((length(L)+1)/2);

%D_k_star(j,i)=L_k_star(j,i)*alpha(i);  %D-T function
D_k_star(j,i)=L_k_star(j,i)*(alpha(i)-lambda(j)+1);  %D-T function
end
end



% % %Calculate D_per
% 
% D_per1=zeros(1,length(alpha_rev));
% D_per2=zeros(1,length(alpha2_vec));
% 
% 
% 
% 
% 
% for i=1:length(lambda)
%     for j=1:length(alpha_rev)
%         D_per1(i,j)=var*(1/(2*alpha_rev(j))-0.5);
%     end
% end
% 
% for i=1:length(var)
%     for j=1:length(alpha2_vec)
%         D_per2(i,j)=(1-alpha2_vec(j))*var(i);
%     end
% end
% 
% 
% D_per=[D_per1 D_per2];

 for j=1:length(lambda)
     plot(alpha,D_k_star(j,:),color(j));
    % plot(alpha,log(D_k_star(j,:)),color(j));
      hold on
%      plot(alpha,D_per(j,:),'r');
 end
 
     set(gca,'PlotBoxAspectRatio',[5 3 1])
     xlabel('$\alpha$','Interpreter','latex');
     ylabel('$D^*_\beta(\alpha)$','Interpreter','latex');



