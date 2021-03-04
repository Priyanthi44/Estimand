
parfor( edge =1:6)
    %     Seletct channel ID
    Channel_ID = selectChannel(edge);
    %initialise parameters
    
    %% RDA
    d=0.99; %negative fb para
    g=0.1;%negative fb para
    n=0; %initial negfb
    p=0;%initial position
    L=10;
    B=2;
    l_max=L; % max proof reading distance
    beta=B; % min threshold
    exceeded =false;
    negfb=false;
    %% TD(lambda)
    
    alpha=0.1;% SGD learning rate
    elig=0; % Eligibility trace
    
    %% Input
    items=10000;
    side =1;
    
    
    
    %% CREATE DATASET
    %
    U=(Dataset(side,2,items));
    
    %  if(~timeseries)
    %      fileID = fopen(strcat('dataset.txt'),'w');
    %         fprintf(fileID,'%12.8f\n',U);
    %         fclose(fileID);
    %  end
    if(~timeseries)
        A =fopen('dataset.txt','r');
        U = fscanf(A,'%f');
    end
    
    
    %% SNN neuron
    dt=0.1;
    t_end= items/10;
    t_vect =0:dt:t_end;
    
    V_vect=zeros(1,length(t_vect));
    I_max=1;
    I_e_vect =ones(1,length(t_vect));
    beta_vect =ones(1,length(t_vect))*8*0.21159;
    V_vect(1)=0;
    max=0.9;
    j=0;
    o=0;
    T=0;
    countsea=0;
    seakeeping=false;
    
    
    %% DNN
    
    inputlayersize=1;
    outputlayersize=1;
    hiddenlayersize=7;
    syn0 = abs(2*rand(inputlayersize,hiddenlayersize) - 1);
    syn1 = abs(2*rand(hiddenlayersize,outputlayersize) - 1);
    k=0;
    %% The Simulation
    
    for i=1:items
        if(~mod(i,15))
            %         Analyse and upload values to the server
            
        end
        if (abs(p)>B)
            exceeded=true;
            if(abs(p)>L)
                seakeeping=true;
            else
                seakeeping=false;
            end
        end
        if(~exceeded)% i.e. p<beta
            if (p>0)
                p= p+U(i)-n;
            else
                p= p+U(i)+n;
            end
            
            n=d*n;
            % Calculate the position
            [V_vect(i+1),I_e_vect(i),j]= runSNN(abs(p),l_max,V_vect(i),j);
        end
        if(seakeeping)
            x=L/B;
            beta=L;
            l_max=abs(p)*x;
        else
            beta=B;
            l_max=L;
        end
        if (exceeded)
            %         %Run DNN
            t=cputime;
            o=o+1;
            
            %% WITH DNN
            [n,syn0,syn1,t,k]= runDNN(n,getOutput(abs(p),U(i),beta,max),syn0,syn1,k,t);
            T=t+T;
            
            %%       WITHOUT DNN
            %        n=d*n+g;
            %        t=cputime-t;
            %        T=t+T;
            %%
            
            if (p>0)
                p= p+U(i)-n;
            else
                p= p+U(i)+n;
            end % Calculate the position
            
            
            %%
            [V_vect(i+1),I_e_vect(i),j]= runSNN(abs(p),l_max,V_vect(i),j);
            if abs(p)<B
                exceeded=false;
            end
            %
        end
        
        
        
        
        if(beta==B)
            countsea=countsea+1;
        end
        
        
    end
    
    %   figure(2);
    %  plot (t_vect,U);
    
    figure(edge);
    
    if(k==0)
        k=o;
    end
    
    T=T/k
    plot (t_vect,V_vect);
    hold on;
    % plot(t_vect,beta_vect, '-r');
    title ('Voltage vs Time');
    xlabel('Time in ms')
    ylabel('Voltagein mV');
end
% Select Channel
function channelID =selectChannel(edge)
channelID_1 = 1;
channelID_2 = 2;
channelID_3 = 3;
channelID_4 = 4;
channelID_5 = 5;
channelID_6 = 6;

    switch (edge)
        case 1
            channelID = channelID_1;        
        case 2
            channelID = channelID_2;
        case 3
            channelID = channelID_3;

        case 4
            channelID = channelID_4;
        case 5

            channelID = channelID_5;

        case 6
            channelID = channelID_6;

    end
end

%% Sensor Data


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


err=abs((p/l)*100); % Percentage of error
I_vect=I_max*(1-err/100); % Adjust input current according to the error
if(err>100)
    I_vect=0;
end
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

if(V_vect>=8*0.21159)
    j=j+1;
    
    
end
if(V_vect<0)
    p
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
function [layer2,syn0,syn1,t,k]=runDNN(trainingdata_inputs,trainingdata_outputs,syn0,syn1,k,t)
k=k+1;
%    error_tolerance = 0.1;
%   if (fail)
%  while (fail)

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
error_tolerance = 0.1;
if errorval<error_tolerance
    fail=false;
end
%print out debug data
%        if iter==1 || mod(iter,100000) == 0
%                 fprintf("\titer=%.0f, Error: %f\n", iter, errorval)
%                 %syn0
%                 %syn1
%        end

%  end
%   else
%
%   end
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
t=cputime-t;
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

