 % This is handled in the backgroung /parallel
% this generated r_loss
function output=runDNN(trainingdata_inputs,inputlayersize,outputlayersize,hiddenlayersize,backpropagate)
persistent syn0
if isempty (syn0)
    syn0=abs(2*rand(inputlayersize,hiddenlayersize) - 1);
    
end
persistent syn1

if isempty(syn1)
    syn1=abs(2*rand(hiddenlayersize,outputlayersize) - 1);
    
end
persistent layer0
if isempty(layer0)
    layer0=0;
    
end
persistent layer1
if isempty(layer1)
    layer1=0;
    
end
persistent layer2
if isempty(layer2)
    layer2=0;
    
end
if(~backpropagate)
  
    layer0=trainingdata_inputs;
    
    layer1=(1)./(1+exp(-1.*(layer0*syn0)));
    %multiply inputs by weights and apply sigmoid activation functoin
    
    layer2=(1)./(1+exp(-1.*(layer1*syn1)));
    output=layer2;
else
    %multiply hidden layer by 2nd set of weights and apply sigmoid activation function
  
    layer2_error=layer2-trainingdata_inputs;
    
    %which direction is the target value
    layer2_delta = layer2_error.*(exp(layer2)./(exp(layer2)+1).^2);
    
    %how much did each l1 value contribute to l2 error
    layer1_error = layer2_delta*syn1.';
    
    %which direction is target l1
    layer1_delta = layer1_error.*(exp(layer1)./(exp(layer1)+1).^2);
    alpha=0.1;
    %adjust values
    %             errorval = mean(abs(layer2_error));
    syn1 = syn1 - alpha.*(layer1.'*layer2_delta);
    syn0 = syn0 - alpha.*(layer0.'*layer1_delta);
    output=0;
end

end