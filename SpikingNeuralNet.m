classdef SpikingNeuralNet < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        noNodes;
        net;
        forces;
        Value;
    end
    
    methods
        function thisNet=SpikingNeuralNet(no)
          if nargin~=0
               
          
            thisNet.noNodes =no;
            thisNet.forces =zeros(1,thisNet.noNodes);
            thisNet.net = Neuron.empty (thisNet.noNodes,0);
%              i=1;
%                 while i~=no+1
%                     % #Encoding the opposite action 
%                    
%                     i=i+1;
%                     
%                 end

         end
               
                                                             
        end 
        function createNet(thisNeuron, i, thisSNN)
            thisSNN.net(i) = thisNeuron;
            
        end
        
        function net= getNet(thisSNN);
            net =thisSNN.net();
        end
        
        function forces= storeForce(force, i)
            forces(i)= force;
        end
        
        function thisNeuron =getNeuron(ind)
            thisNeuron  = Neuron.empty ;
            thisNeuron = thisNet.net(ind);
    end
    
end
end
