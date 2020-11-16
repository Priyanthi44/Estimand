function Exp(Attribute)

%initialise parameters
 %% RDA
 d=0.99; %negative fb para
 g=0.1;%negative fb para
 n=0; %initial negfb
 p=0;%initial position
 l_max=10; % max proof reading distance
  beta=2; % min threshold
 exceeded =false;
 negfb=false;
 %% TD(lambda)
 gamma=0.99; % RL discount factor
 lamda=1; % RL lambda
 alpha=0.1;% SGD learning rate
 elig=0; % Eligibility trace
 
 %% Input
 items=5000;
 side =1;

%  for i=1:6
%      if (mod(i,2))
%          if (Attribute(i)>0)
%              neg_axis=Dataset(side,Attribute(i),items);
%          else
%             neg_axis=0; 
%          end
%      else
%          if (Attribute(i)>0)
%              pos_axis=Dataset(side,Attribute(i),items);
%          else
%              pos_axis=0;
%          end
%      end
%      if i==2
%          if neg_axis==0
%              x_axis=-neg_axis;
%          elseif pos_axis ==0
%               x_axis=pos_axis;
%          else
%               x_axis=pos_axis-neg_axis;
%          end
%      elseif i==4
%          if neg_axis==0
%              y_axis=-neg_axis;
%          elseif pos_axis ==0
%               y_axis=pos_axis;
%          else
%               y_axis=pos_axis-neg_axis;
%          end
%      elseif i==6
%          if neg_axis==0
%              z_axis=-neg_axis;
%          elseif pos_axis ==0
%               z_axis=pos_axis;
%          else
%               z_axis=pos_axis-neg_axis;
%          end
%      end
%      
%  end
%  U= sqrt((x_axis).^2+(y_axis).^2+(z_axis).^2);
%  A=U;
%  A(A==0) =1;
%  ang_x=x_axis/U;
%  ang_y=y_axis/U;
%  ang_z=z_axis/U;
 datatype='norm';
 if strcmp(datatype,'const')
     data=1;
 elseif strcmp(datatype,'uni')
         data=2;
 elseif strcmp(datatype,'norm')
         data=3;
 elseif strcmp(datatype,'sin')
         data=4;
 elseif strcmp(datatype,'lin')
         data=5;
 elseif strcmp(datatype,'stp')
         data=6;
 elseif strcmp(datatype,'exp')
         data=7;
 else 
         data=8;
 end
 
%   
U=Dataset(side,data,items)

%% SNN neuron
        dt=0.1;
        t_end= items/10;
        t_vect =0:dt:t_end;
       % t_items=1:1:items;
        V_vect=zeros(1,length(t_vect));
        I_max=1;
        I_e_vect =ones(1,length(t_vect));
        beta_vect =ones(1,length(t_vect))*8*0.21159;
         V_vect(1)=0;
        max=0.9;  
        j=0;
        o=0;
        c=0;
        d=0;
   % figure(2);
      [counts, bins] = hist(U);
   %plot(bins, counts);
   mean(U)
   std(U)%# get a line plot of the histogram     
  %% DNN 
       
        inputlayersize=1;
        outputlayersize=1;
        hiddenlayersize=7;
        syn0 = abs(2*rand(inputlayersize,hiddenlayersize) - 1);
        syn1 = abs(2*rand(hiddenlayersize,outputlayersize) - 1);
        k=0;
  %% The Simulation
       
 for i=1:items          
   if (p>beta)
      exceeded=true; 
%       if(negfb==false)% initialise n only once 
%           n=d*n+g;
%           negfb=true;
%       end
   end
  if(~exceeded)% i.e. p<beta
     p= p+U(i)-n;
     n=d*n;
     % Calculate the position  
    [V_vect(i+1),I_e_vect(i),j]= runSNN(p,l_max,V_vect(i),j);
  end
    
    if exceeded
%         %Run DNN
        
         o=o+1;
          [n,syn0,syn1,k]= runDNN(n,getOutput(p,U(i),beta,max),syn0,syn1,k);
        
% if((size(n)==1)& (size(d)==1) &(size(g)==1))
        p= p+U(i)-n; % Calculate the position 
        [V_vect(i+1),I_e_vect(i),j]= runSNN(p,l_max,V_vect(i),j);
         if p<beta
         exceeded=false;
         end
%       %Run Kalman Filter
%         [priori, posteriori] = KalmanFilter(U(i),(d+Heaveside(p,beta,g,n)),getPosition(I_max,l_max,I_e_vect(i)),n); 
        %Estimate current state and future state after applying negfb
%      
%      r=calculateReward(posteriori,priori);%Calculate future reward
%     reward = r*ones(1,7);
% 
%     %Update theta
% v1=gamma.*syn0.*posteriori+reward;
% v2=syn0.*priori;
%     delta_1= v1-v2 ; 
%     v1=gamma.*syn1.*posteriori+transpose(reward);
%     v2=syn1.*priori;
%     delta_2= v1-v2 ;  
%     
%     
%     elig =gamma.*lamda.*elig+priori;
% 
%     syn0 = syn0+alpha.*elig.*delta_1;
%     syn1 = syn1+alpha.*elig.*delta_2;
% end
    end
    
 
     if(I_e_vect(i)>=0.999)
        % n=0;
         d=d+1;
       
         
     end

 
 end
 
%   figure(2);
%  plot (t_vect,U);

  figure(1);
  j
  k
  o
  d
  c
        plot (t_vect,V_vect);
        hold on;
       plot(t_vect,beta_vect, '-r');
        title ('Voltage vs Time');
        xlabel('Time in ms')
        ylabel('Voltagein mV');
        
        
%% Sensor Data   

%     function U=getSideData(side,data,items)
%      %  datatype='lin';
%  if strcmp(datatype,'const')
%      data=1;
%  elseif strcmp(datatype,'uni')
%          data=2;
%  elseif strcmp(datatype,'norm')
%          data=3;
%  elseif strcmp(datatype,'sin')
%          data=4;
%  elseif strcmp(datatype,'lin')
%          data=5;
%  elseif strcmp(datatype,'stp')
%          data=6;
%  elseif strcmp(datatype,'exp')
%          data=7;
%  else 
%          data=8;
%  end
%  
% %   
% U=Dataset(side,data,items);   
    end
function p=getPosition(I_max,l_max,I)
current=I/I_max;
current=1-current;
        p=l_max*current;
end
%% SNN 
function [V_vect, I_vect,j]= runSNN(p,l,V_vect_1,j)  
        dt=0.1;
        R_m =10;
        tau =100;
        Th=8;
        I_max=1;
        V_th=Th;
        V_reset = 0;
       
        
        err=(p/l)*100; % Percentage of error 
       I_vect=I_max*(1-err/100); % Adjust input current according to the error
       
    V_inf = I_vect* R_m;
    V_vect =V_inf + (V_vect_1-V_inf) *exp(-dt/tau); %Calculate voltage
%    Select threshold voltage to reflect the error percentage
    if (err>0 && err<3.33)
             V_th=Th*0.74869;
   elseif (err>3.33 && err<15.00)
            V_th=Th*0.48828;
   elseif (err>15.00 && err<26.67)
            V_th=Th*0.22786;
   elseif (err>26.67 && err<36.67)
            V_th=Th*0.21159;
   elseif (err>36.67 && err<46.67)
            V_th=Th*0.19531;
   elseif (err>46.67 && err<53.33)
            V_th=Th*0.17903;
   elseif (err>53.33 && err<100.00)
            V_th=Th*0.16276;
   end    

   if(V_vect>=Th)
         j=j+1;
         
   end
 
   %indicate spike
                if (V_vect>V_th)
                      
                    V_vect =V_reset;
                    
                end

         
  end
        
%% Calculate Reward        
 function [n]=getOutput(p,u,beta,max)
 n= p+u-beta;
 if n>max
   n=max-u;  
 end

 end       
function r=calculateReward(p,pr)
if(p<pr)
    r=1;
elseif(p>pr)
    r=-1;
else
    r=0;
end
end
%% Run DNN
function [layer2,syn0,syn1,k]=runDNN(trainingdata_inputs,trainingdata_outputs,syn0,syn1,k)
   k=k+1;
    layer0=trainingdata_inputs;
    
    layer1=(1)./(1+exp(-1.*(layer0*syn0))); 
    %multiply inputs by weights and apply sigmoid activation functoin
    
    layer2=(1)./(1+exp(-1.*(layer1*syn1)));
    %multiply hidden layer by 2nd set of weights and apply sigmoid activation function
    
    %cost function (how much did we miss)
    layer2_error=layer2-trainingdata_outputs;
    
    %which direction is the target value
    layer2_delta = layer2_error.*(exp(layer2)./(exp(layer2)+1).^2);

        %how much did each l1 value contribute to l2 error
     layer1_error = layer2_delta*syn1.';

        %which direction is target l1
     layer1_delta = layer1_error.*(exp(layer1)./(exp(layer1)+1).^2);
     alpha=0.1;
     %adjust values
            errorval = mean(abs(layer2_error));
            syn1 = syn1 - alpha.*(layer1.'*layer2_delta);
            syn0 = syn0 - alpha.*(layer0.'*layer1_delta);
            
            %print out debug data
%        if iter==1 || mod(iter,100000) == 0
%                 fprintf("\titer=%.0f, Error: %f\n", iter, errorval)
%                 %syn0
%                 %syn1
%        end
       error_tolerance = 0.05;

%         if errorval<error_tolerance
%         
%              fprintf("Stopping at: %f error\n", errorval)
%               
%         
%         end
%         layer0 =n;
% layer1=max(0,(layer0*syn0+syn1));%applying reLu
% %  layer2=min(1,layer1);%filtering :applying constraint (max distance)
% % layer2(layer2==max(layer2(:)))=0;
% out=max(layer1);
% if(size(out)==1)
% index=find(layer1==out);
% d=syn0(index);
% g=syn1(index);
% else
%     out=n;
% end
end
        
%% Get parameter g        
function g=Heaveside(p,beta,g,n)
if(p>beta)
    if(abs(n)>0)
        g=g/n;
    else
        g=0;
    end
            
else
    g=0;
end
        
            
end
%% Run Kalman Filter        
function [priori,posteriori]= KalmanFilter(u,para,p,n)

A = [1 -1; 0 para] ; % state transition matrix:  
B = [1; 0]; %input control matrix
C = [1 0]; % measurement matrix

Q= [p; n]; %initized state--it has two components: [position; negative feedback] 
Q_estimate = Q;  %x_estimate of initial location 
position_noise_mag = 0.05; %process noise: the variability in position measurements
negativefb_noise_mag = 0.3;  %measurement noise:stdv of location, in meters
Ez = position_noise_mag^2;% Ez convert the measurement noise (stdv) into covariance matrix
Ex = [position_noise_mag^2 position_noise_mag*negativefb_noise_mag; position_noise_mag*negativefb_noise_mag negativefb_noise_mag^2]; % Ex convert the process noise (stdv) into covariance matrix
P = Ex; 
% Predict next state of the quail with the last state and predicted motion.
    Q_estimate = A * Q_estimate + B * u;
     %predic_state =  Q_estimate(1) ;
     priori=  Q_estimate(1) ;
    %predict next covariance
    P = A * P * A' + Ex;
%      predic_var =  sqrt(P(1)) ;
%     % predicted Ninja measurement covariance
%     range=  p-5:.01:p+5;
%     priori= normpdf(range,predic_state,predic_var);

    % Kalman Gain
    K = P*C'*inv(C*P*C'+Ez);
    % Update the state estimate.
    Q_estimate = Q_estimate + K * (p - C * Q_estimate);
    % update covariance estimation.
    %P =  (eye(2)-K*C)*P;
    %Store for plotting
     %posteriori= normpdf(range,Q_estimate(1),P(1));
     posteriori=Q_estimate(1);
%     Q_loc_estimate = [Q_loc_estimate; Q_estimate(1)];
%     vel_estimate = [vel_estimate; Q_estimate(2)];
%     P_mag_estimate = [P_mag_estimate; P(1)]
end
 