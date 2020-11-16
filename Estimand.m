% Run the DRL system
%
%If the postion is within safe limits -> RDA-> SNN
%if the position is not within safe limits ->DNN->SNN
%
% 
%The mathematical model-Feeds the DNN response from a seperate thread
% The mathematical model access the current response <- from eq/ system
% constraints->p,e,beta/ weights syn0 syn1/


% Global /Parallel parameters
 % weights
 %
%use eq to calculate n_t
%calculate n_t+1 using 
%forward propogate DNN ->get n_t+1
%backpropagete DNN 
%calculate loss
%optimise parameters with n_t and loss
