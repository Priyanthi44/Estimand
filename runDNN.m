% This is handled in the backgroung /parallel
% this generated r_loss
function r_t=runDNN()
layer0=trainingdata_inputs;

layer1=(1)./(1+exp(-1.*(layer0*syn0)));
%multiply inputs by weights and apply sigmoid activation functoin

r_t=(1)./(1+exp(-1.*(layer1*syn1)));
%multiply hidden layer by 2nd set of weights and apply sigmoid activation function

end