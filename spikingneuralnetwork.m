function spikingneuralnetwork(items,Attributes,algo,impact,adj)
%% Attributes consists of 6x3 matrix of force type, intensity,vertical or not, for each side 
%% TODO items*1000

 % one node for each size
noNodes= 6; 
noofitems=items*1000; % iterations
% minimum force required to move the vessel
%% Constant values such as minimum forces and limits for both horizontal and vertical dimentions
%%
 t =TCR.const;
% b=t.beta;
% m=t.mass;
% l1= t.limit1;
% l2=t.limit2;
% l=t.limit;
%% create the spiking neural network with six nodes
snn = SpikingNeuralNet(noNodes);

            
for j=1:noNodes
     if mod(j,2)
          
          thisNeuron=Neuron(j+1,items);
          createNet(thisNeuron,j,snn);
     else
         thisNeuron=Neuron(j-1,items);
          createNet(thisNeuron,j,snn);
     end
    
end
%% run neurons
result=cell(6,1);
voltage=zeros(noofitems,1);

for k=1:noNodes
    result{k}=runNeuron(k,snn.net(k),Attributes(k),Attributes(k+6),Attributes(k+12),snn,impact,adj);
%     voltage=snn.net(k).V_vect +voltage;
%     result{k}(result{k}<0)=0;
end    
%  figure(1);
%         plot (TCR.const.time,voltage);
%         title ('Voltage vs Time');
%         xlabel('Time in ms')
%         ylabel('Voltagein mV');  

%% parallally run each node
% parallel.defaultClusterProfile('local');
% c=parcluster();
% obj = createJob(c);
% createTask(obj,@runNeurons,1,{{snn.net(1),Attributes(1),Attributes(7),Attributes(13),snn} {snn.net(2),Attributes(2),Attributes(8),Attributes(14),snn} 
%     {snn.net(3),Attributes(3),Attributes(9),Attributes(15),snn} {snn.net(4),Attributes(4),Attributes(10),Attributes(16),snn} 
%     {snn.net(5),Attributes(5),Attributes(11),Attributes(17),snn} {snn.net(6),Attributes(6),Attributes(12),Attributes(18),snn}});
% submit(obj);
% wait(obj);
% result= fetchOutputs(obj);

%%
%% Get the resultant force applied

if(algo==1)
F=zeros(noofitems.*2,noNodes);
P=zeros(2*noofitems,noNodes/2);
S=zeros(2*noofitems,noNodes/2);
else
F=zeros(noofitems,noNodes);
P=zeros(noofitems,noNodes/2);
S=zeros(noofitems,noNodes/2);
end

for i=1:noNodes
    if exist(strcat( num2str(i),'.txt'), 'file')
        fileID = fopen(strcat( num2str(i),'.txt'),'r');
        formatSpec = '%f';
        A = fscanf(fileID,formatSpec);
        fclose(fileID);
        %delete(strcat( num2str(i),'.txt'),'file');
% multiply with attribute 
            vertical =Attributes(i+12);
            intensity=Attributes(i+6);
            
            if(impact==1)
            maxforce= TCR.getMax(vertical);
            else
            maxforce= TCR.getMin();  
            end
            minforce=plus(intensity,1);
            minforce =times(minforce,maxforce);
               % 
% selects conditions that the vessel is moving
        if i<5
            for j=1:size(A,1)
                A(j)=minforce .*A(j);
                
%             A(j)=minus(min,b);
            
            end
        end
    A(A<0) =0;
  z=zeros(noofitems.*2,1);
  k=1;
  for j=1:2*noofitems
      if(mod(j,2))
      z(j)=A(k);
      else
          z(j)=-result{i}(k);
          k=k+1;
      end
  end
  
size(A)
size(result{1})
%  factor =intensity/intensity^2;
if(algo==2)
 F(:,i) = A;
%  result{i} =(factor).*result{i};
else
%   F(:,i) = A-result{i};
F(:,i) = z;
  
        

end

    end
end

j=1;

for i=1:noNodes
     if mod(i,2)% odd number
           P(:,j)=(F(:,i)-F(:,i+1)); %calculate force of the whole column
           j=j+1;
     end
end 

tic
for j=1:noNodes/2
    
   for i=1:noofitems
        
        if(i==1)
        v=P(i,j).*(1/1000)/TCR.const.mass;
        else
            v =u+P(i,j).*(1/1000)/TCR.const.mass;
        end
        
        S(i,j)= getDistance(P(i,j),v);
        u=v;
%         if(~mod(i,100))
%         
%         else
%            if(i>1)
%            S(i,j)= S(i-1,j);  
%            else
%               S(i,j)= 0; 
%            end
%         end
%        if(i>1)
%           
%       S(i,j)= (2*P(i-1,j)- P(i,j))/(2*m*1000*1000);
%        else
%           S(i,j)= P(i,j)/(2*m*1000*1000);
%        end
   end
end
toc
t=toc;
    
%  S(S<0) =0;
Simulate(S);
% destroy(obj);
% clear
% delete(gcp)
% for i=1:size(snn.net,2)
%   
% 
% end
%  n=Neuron(i,i+1);
% force=runNeuron(n,i);
% nn.forces(i) =force;
    
% Responses(NeuronMatrix,i,i,i+1);

end
function distance=getDistance(f,u)
        %s=ut+1/2at2
        distance= u.*(1/1000)+1/2.*(f/TCR.const.mass).*(1/1000)^2;
        
 end
% initiateNeuron(NeuronMatrix(2));
% runNeuron(NeuronMatrix(2),2);
% parfor i=1:num
%     spikingneuron(i)
% end

 
%   end

% end
% num=2;
function  thrust=runNeurons(side,node, n,intensity,vertical, net,impact)

   
    thrust=runNeuron(side,node,n,intensity,vertical,net,impact);
   
 
   

end
% force=net.force;
% thrust=node.thrust ;



function nArray= getArray(i,SingleLayer)
count=0;
    
 for j=1:(size(SingleLayer,2))
   
        if(SingleLayer(j)~=SingleLayer(i))
            count=count+1;
        nArray(count)=SingleLayer(j);
        end
        
           
 end
 
end
