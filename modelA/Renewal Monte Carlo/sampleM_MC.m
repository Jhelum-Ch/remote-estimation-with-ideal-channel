function [count,T, D, iter] = sampleM_MC(p,k,epsilon,ITERW_space = [-1, 0, 1];
W_PMF   = [p, 1-2*p, p];

X = 0; t=0;
D = 0; T = 1;
%scale_beta=1;
iter=1;
count=1;
%while (abs(X) < k || rand(1) < epsilon) && iter < ITER
while abs(X) < k && iter < ITER
%while iter < ITER
  W = randsample(W_space, 1, true, W_PMF);
  %scale_beta=scale_beta*beta;
  t=t+1;
    Xnext = X + W;
  if (abs(Xnext) >= k)
      %count=count+scale_beta;
      return;
  end
  T = T + 1;
  D = D + abs(Xnext);
%    if (abs(X) >= k && rand(1) > epsilon)
%        Xnext = W;
%    else
%        Xnext = X + W;
%    end

  X = Xnext;
  iter = iter+1;
end