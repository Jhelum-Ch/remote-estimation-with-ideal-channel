clear all
clc
format long


lambda=1;




% %Code for fixed K, different variances%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% K=10;
% 
% a=-K; b=K;
% 
% p=0.3;
% 
% 
% 
% 
% Var=0.5:0.8:5;
% nfinal=zeros(length(Var),1);
% Color=['b';'r';'k';'g';'m';'c'];
% 
% for i=1:length(Var)
% var=Var(i);
% 
% kernel=@(e,x) (1/sqrt(2*pi*var))*exp(-(x-e).^2/(2*var));
% 
% 
% RHS=@(e) e.^0;
% 
% AbsTol = 1e-6;
% RelTol = 1e-3;
% 
% [sol,errest,cond] = Fie(lambda,a,b,1,kernel,RHS,AbsTol,RelTol);
% M=sol.x;
% t = sol.s;
% nfinal(i) = length(t) - 1;
% 
% % Interpolate the solution for assessing the error and if necessary,
% % use to get a smooth graph.
% tint = linspace(a,b,150);
% xint = ntrpFie(sol,tint);
% 
% 
% if nfinal(i) < 150
%     plot(tint,xint,Color(i));
%     set(gca,'PlotBoxAspectRatio',[5 3 1])
%     xlabel('$e$','Interpreter','latex');
%     ylabel(' $L^{(k)}(e)$','Interpreter','latex');
%     hold on
% else
%     plot(t,M,Color(i));set(gca,'PlotBoxAspectRatio',[5 3 1])
%     xlabel('$e$','Interpreter','latex');
%     ylabel(' $L^{(k)}(e)$','Interpreter','latex');
%     
%     hold on
% end
% end














%Code for a known set of values of K

K=linspace(0.1,100);
a=zeros(length(K),1);
b=zeros(length(K),1);

for i=1:length(K)
a(i)=-K(i); b(i)=K(i);
end

p=0.3;

%var=2*p;
%var=10;
var=[0.5;1.0;10];

color=['b';'r';'k'];




RHS=@(e) e.^0;
RHS1=@(e) e.^2;



AbsTol = 1e-6;
RelTol = 1e-3;
M_k=zeros(length(var),length(K));
L_k=zeros(length(var),length(K));
D_k=zeros(length(var),length(K));
N_k=zeros(length(var),length(K));


for j=1:length(var)
    kernel=@(e,x) (1/sqrt(2*pi*var(j)))*exp(-(x-e).^2/(2*var(j)));
    
for i=1:length(K)
    M=[];
    L=[];
    t1=[];
    t2=[];
[sol1,errest1,cond1] = Fie(lambda,a(i),b(i),1,kernel,RHS,AbsTol,RelTol);
[sol2,errest2,cond2] = Fie(lambda,a(i),b(i),1,kernel,RHS1,AbsTol,RelTol);
M=sol1.x;
t1 = sol1.s;
M_k(j,i)=M((length(M)+1)/2);
L=sol2.x;
t2 = sol1.s;
L_k(j,i)=L((length(L)+1)/2);

D_k(j,i)=L_k(j,i)/M_k(j,i);
N_k(j,i)=1/M_k(j,i);




nfinal1 = length(t1) - 1;

% Interpolate the solution for assessing the error and if necessary,
% use to get a smooth graph.
tint1 = linspace(a(i),b(i),150);
xint1 = ntrpFie(sol1,tint1);


nfinal2 = length(t2) - 1;

% Interpolate the solution for assessing the error and if necessary,
% use to get a smooth graph.
tint2 = linspace(a(i),b(i),150);
xint2 = ntrpFie(sol2,tint2);


end


end
% 


figure

for j=1:length(var)
    
plot(K,M_k(j,:),color(j));
set(gca,'PlotBoxAspectRatio',[5 3 1])
xlabel('$k$','Interpreter','latex');
ylabel(' $M^{(k)}(0)$','Interpreter','latex');
hold on
end


figure
for j=1:length(var)


plot(K,L_k(j,:),color(j));
set(gca,'PlotBoxAspectRatio',[5 3 1])
xlabel('$k$','Interpreter','latex');
ylabel(' $L^{(k)}(0)$','Interpreter','latex');
hold on
end


figure
for j=1:length(var)

plot(K,D_k(j,:),color(j));
set(gca,'PlotBoxAspectRatio',[5 3 1])
xlabel('$k$','Interpreter','latex');
ylabel(' $D^{(k)}(0)$','Interpreter','latex');
hold on
end


figure
for j=1:length(var)

plot(K,N_k(j,:),color(j));
set(gca,'PlotBoxAspectRatio',[5 3 1])
xlabel('$k$','Interpreter','latex');
ylabel(' $N^{(k)}(0)$','Interpreter','latex');
hold on
end
