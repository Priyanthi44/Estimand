rrrreeeeeeeeeeeeerrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrz\wclassdef Neuron <SpikingNeuralNet
    %NEURON Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        items;
        I_e_vect;
        opposite;
         thrust;
         force;
         inputforce;
         t_vect;
         V_vect;
    end
    
    properties (Access=protected)
        dt=0.001;
        R_m =10;
        tau =100;
        V_th=8;
        V_spike =10;
        V_reset = -1;
        t_end;
        
       
        
    end
    
    properties (Constant) 
        Th=8;
    end 
    
    properties (Dependent)
    end
    
    methods(Static)
          function x=setItems()
            x= 500;
        end 
        function x =setEnd(y)
            x= y/10;
        end 
    end 
    
    methods
        function thisNeuron=Neuron(opp,noofitems)
            thisNeuron = thisNeuron@SpikingNeuralNet();
                
            if nargin~=0
                
                thisNeuron.opposite =opp;
                thisNeuron.items= noofitems-thisNeuron.dt;
                thisNeuron.t_end =thisNeuron.items;
                thisNeuron.t_vect =0:thisNeuron.dt:thisNeuron.t_end;
                t=TCR.const;
                t.time=thisNeuron.t_vect;
                thisNeuron.V_vect=zeros(1,length(thisNeuron.t_vect));
                thisNeuron.thrust=0;
                thisNeuron.force=0;
                thisNeuron.inputforce=0;
                thisNeuron.I_e_vect = ones(1,length(thisNeuron.t_vect));
                thisNeuron.items= size(thisNeuron.t_vect,2);
                thisNeuron.V_vect(1)=0;
             
              
            end 
        end
function inihibitoryResponse(thisNet, response,opposite, time,items)
    % find the opposite neuron
    
     oppoNeuron =thisNet(opposite);
     %check if opposite neuron has a force
    % check if I is not equal to one
    if (oppoNeuron.I_e_vect(time)<1)
        % if so increase the current
        if(oppoNeuron.I_e_vect(time)+response<1)
            if(items<(time+1))
        oppoNeuron.I_e_vect(time+1) = oppoNeuron.I_e_vect(time)+response;
        else
          oppoNeuron.I_e_vect(time+1) = 1; 
            end
        end
    end%# reflecting the force
%     end
    
         
           
end
     
function  thrustresponse = runNeuron(side,thisNeuron,type,intensity,vertical,snn,impact,adj) % k is the choice of dataset
            r= TCR(vertical); % create a T-cell receptor
           
            spikes=0;
           
            if(impact==1)
            maxforce= TCR.getMax(vertical);
            else
            maxforce= TCR.getMin();  
            end
            minforce=plus(intensity,1);
            minforce =times(minforce,maxforce);
            
            Y=Dataset(side,type,thisNeuron.items);
            j=0;
            
            thrustresponse =zeros(1,length(thisNeuron.t_vect));
            for i=1:thisNeuron.items
                
                u=Y(i);% input environmental force
                
                
                u=u*minforce;
                thisNeuron.force=thisNeuron.force+u;
                j=j+1;
                if (u<0)
                    
                    u=0;
                end
%                 thisNeuron.force(i)= u; 
                
                y=1;
                 if (thisNeuron.V_th~=Neuron.Th)% # Spiking has been reduced
                     
                     y=1/(Neuron.Th-thisNeuron.V_th);
                     
                     
%                      u=u-thisNeuron.force(i-1);% removing the force from the thrusters
                     %call other neurons
                    inihibitoryResponse(snn.net,y,thisNeuron.opposite,i,thisNeuron.items);
%                     thisNeuron.thrust(i)= thisNeuron.response;   
                    u=u-thisNeuron.thrust;
                    if (u<0)
                        u=0;
                    end
                 
                 else
                     
                 end 
                 
                 x=runRDA(r,u,vertical); % run RDA 
                
                
                if(x>0 && x<1)
                    thisNeuron.inputforce=thisNeuron.force-(j.*maxforce);
                   if(adj==1)
                 thisNeuron.thrust=  (x+1).*maxforce; 
                    else
                  thisNeuron.thrust=  (maxforce+thisNeuron.inputforce).*(1+x+intensity); 
                  thisNeuron.inputforce=0;
                  thisNeuron.force=0;
                   j=0;
                    end
                   thisNeuron.I_e_vect(i) =(thisNeuron.I_e_vect(i)-x);
                   
                  
      % tuning the firing threshold according to the disturbance from RDA 
                     if x>0.75 && x<1
                        thisNeuron.V_th = Neuron.Th*0.125;
                    elseif x>0.5 && x<0.75
                        thisNeuron.V_th = Neuron.Th*0.375;
                    elseif x>0.25  && x<0.5
                        thisNeuron.V_th = Neuron.Th*0.625;
                    elseif  x<0.25 
                        thisNeuron.V_th = Neuron.Th*0.875;
                    end 
         
       % change the threshold according to the current reduction
                else

                end
                
%   
                V_inf =thisNeuron.I_e_vect(i) * thisNeuron.R_m;
                thisNeuron.V_vect(i+1) =V_inf + (thisNeuron.V_vect(i)-V_inf) *exp(-thisNeuron.dt/thisNeuron.tau);
  
                if (thisNeuron.V_vect(i+1)>thisNeuron.V_th)

                    thisNeuron.V_vect(i+1) =thisNeuron.V_reset;
                    
                    spikes=spikes+1;

               
                end
 
              thrustresponse(i)=thisNeuron.thrust;
             
 
            end
            
% 
         %clear Responses.exhibitoryResponse;
% %     drawnow;
          
         
    
          
      end 
        
    function t_isi=initiateNeuron(thisNeuron)
        spikes=0;
        
        for i=1:thisNeuron.items
             V_inf =thisNeuron.I_e_vect(i) * thisNeuron.R_m;
             thisNeuron.V_vect(i+1) =V_inf + (thisNeuron.V_vect(i)-V_inf) *exp(-thisNeuron.dt/thisNeuron.tau);
  
                if (thisNeuron.V_vect(i+1)>thisNeuron.V_th)
                    thisNeuron.V_vect(i+1) =thisNeuron.V_reset;
                    thisNeuron.V_plot(i+1) =thisNeuron.V_spike;
                    spikes=spikes+1;
                else
                    thisNeuron.V_plot(i+1)=thisNeuron.V_vect(i+1);
                end
        end

        t_isi =5000/spikes;
        drawnow;
        figure(1);
        plot (thisNeuron.t_vect,thisNeuron.V_plot);
        title ('Voltage vs Time');
        xlabel('Time in ms')
        ylabel('Voltagein mV');  
        hold on;
    end
     
 
    function displayNeuron(thisNeuron)
%          thisNeuron.R_m =10;
%                 thisNeuron.tau 
%                 thisNeuron.dt
%                 thisNeuron.V_th
%                 thisNeuron.V_spike 
%                 thisNeuron.V_reset 
%                 thisNeuron.items
%                 thisNeuron.t_end 
%                 
% %                 
% %                 thisNeuron.t_vect
% %                 
% %                 thisNeuron.V_vect
% %                 thisNeuron.V_plot 
%                 
%                 
                thisNeuron.I_e_vect 
%          
%                 thisNeuron.V_vect(1)
%                 thisNeuron.V_plot(1)
    end
        
    end
end

