# remote-estimation-with-ideal-channel
This repository contains MATLAB codes for remote estimation with two agents connected via an idea channel. We consider two models, Model A and Model B. The codes for each model is kept in a separate folder.

Model A: A symmetric scalar Birth-Death Markov chain with countable state space; the birth rate is denoted by 'p'. From any state 'x', the transition probabilities to states 'x+1', 'x-1' and 'x' are 'p', 'p' and '1-2p' respectively. 

Model B: First-order scalar autoregressive process with zer-mean Gaussian innovations. 

These codes reproduce the results presented in the paper:

J. Chakravorty and A. Mahajan, "Fundamental limits of remote estimation of autoregressive Markov processes under communication constraints", IEEE Transactions on Automatic Control, March 2017.
