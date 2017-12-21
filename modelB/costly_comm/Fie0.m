function R = Fie0(beta,k,kernel,RHS1,AbsTol,RelTol)

rescaled_RHS1=@(x) RHS1(x)/beta;
scalar=1/beta;

    %sol = Fie(scalar,-k,k,1,kernel,rescaled_RHS1,AbsTol,RelTol);
    sol = Fie(scalar,-k,k,1,kernel,rescaled_RHS1,AbsTol,RelTol);
    Rvec = sol.x;
    R = Rvec((length(Rvec)+1)/2); 
end