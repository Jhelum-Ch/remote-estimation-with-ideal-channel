clear all
clc

format long

p=0.3;

beta=[0.9;0.95;1.0];

c=0.1;
c_final=20;
c_beta=c;
lambda=zeros(length(beta),1);
 
for i=1:length(beta)
    D=-2-(1-beta(i))/(beta(i)*p);
         lambda(i)=acosh(-D/2);

        k=1;
        c=0.1;
    if beta(i)==1
        c_beta=[];
        L=[];
        M=[];
     while k>0 && c<c_final  
     c_beta(k)=(k*(k+1)*(k^2+k+1))/(6*p*(2*k+1));
     [L(i,k),M(i,k)]=calcLM(k,lambda(i),beta(i),p);
     c=c_beta(end);
     J(k) = (L(i,k)+c_beta(k))/M(i,k);
     k=k+1;
     end
     c_beta1=c_beta;
    elseif beta(i)==0.9
        c_beta=[];
        L=[];
        M=[];
        

         
         while k>0 && c<c_final

         [L(i,k),M(i,k)]=calcLM(k,lambda(i),beta(i),p);
         D(i,k)=L(i,k)/M(i,k);
         if k>1 && c<c_final
         c_beta(k-1)=(D(i,k) - D(i,k-1))/(1/M(i,k-1) - 1/M(i,k));
         c=c_beta(k-1);
         J_disc09(k-1) = L(i,k-1)/M(i,k-1)+ c*((1/M(i,k-1))-(1-beta(i)));
         c_beta09(k-1)=c_beta(k-1); 
         end
            k=k+1;
         end
    else  c_beta=[];
        L=[];
        M=[];
        

         
         while k>0 && c<c_final
             [L(i,k),M(i,k)]=calcLM(k,lambda(i),beta(i),p);
             D(i,k)=L(i,k)/M(i,k);
         if k>1 && c<c_final
         c_beta(k-1)=(D(i,k) - D(i,k-1))/(1/M(i,k-1) - 1/M(i,k));
         c=c_beta(k-1);
         J_disc095(k-1) = L(i,k-1)/M(i,k-1)+ c*((1/M(i,k-1))-(1-beta(i)));
         c_beta095(k-1)=c_beta(k-1); 
         end       
            k=k+1;
         end
    end


end



%Plots
C=[1, 5, 10, 15, 20, 25, 30, 35];

plot(c_beta09,J_disc09)
hold on
plot(c_beta095,J_disc095,'r')
hold on
plot(c_beta1,J,'k')

set(gca, 'XTick', C); % Change x-axis ticks
set(gca, 'XTickLabel', C); 
set(gca,'PlotBoxAspectRatio',[5 3 1])
xlabel('$\lambda$','Interpreter','latex');
ylabel(' $C^*_\beta(\lambda)$','Interpreter','latex');

