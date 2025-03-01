function [irf,gradIrf] = genIRFShockRank(q,vma,neg)
% Generate horizon-h impulse response of the ith variable to the first 
% structural. Also generates analytical gradient.
% Inputs:
% q: first column of Q
% vma: row vector of coefficients in orthogonal reduced-form VMA 
% representation for ith variable at hth horizon
% n: number of variables
% neg: return negative (positive) impulse response if neg = -1 (1)

irf = neg*vma*q; % IRF of ith variable to first shock
gradIrf = neg*vma'; % Gradient 

end