function [L,M]=calcLM(k,m_beta,beta,p)
    a1=sinh(k*m_beta);
    a2=sinh(m_beta);
    a3=sinh(m_beta/2);
    a4=sinh(m_beta);
    a5=cosh(k*m_beta);
    a6=sinh(k*m_beta/2);
    
    if beta==1
        L=(k*(k^2-1))/(6*p);
        M=(k^2)/(2*p);
    else
    L=(a1-k*a2)/(4*beta*p*(a3)^2*a4*a5);
    M=((a6)^2)/(2*beta*p*(a3)^2*a5);
    end

