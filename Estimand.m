% Run the DRL system
%spmd single program multiple data spmd can be used to share data between
%functions
% ****Use spmd if you require communication between workers during a
% computation.*** We do not require communication between workers during a
% computation
%Because we only access the output of each function to process
% labindex- unique index (used in switch/if else)of a given worker numlabs-number of workers
% labsend to specific worker @param (data, labindex))
%labreceive to receive data from(labindex)-labsend and labreceive has to be
%coupled
%
%****Use parfeval or parfevalOnAll if your code can be split into a set of tasks,
%where each task can depend on the output of other tasks.*****
% Data needs to be shared between runDNN and optimisePara and backPropagate
%
% Data needs to send from the runDNN is r_t
%
% spmd block has to be put in
% parpool(2);
% f1 = parfeval(@Function_I, 1);
% f2 = parfeval(@Function_II, 1);
% % fetchNext waits for one of the functions to complete,
% % and also gets the result
% [idx, result] = fetchNext([f1, f2]);
% % We're done, so we can cancel f1 and f2 (one will actually already be complete)
% cancel([f1, f2]);
% parpool(2);

% spmd
%   done = false;
%   state = []; % state used by Function_I or Function_II
%   while ~done
%     % Run Function_I/Function_II for a while
%     if labindex == 1
%       [state, gotSolution] = Function_I(state);
%     elseif labindex == 2
%       [state, gotSolution] = Function_II(state);
%     end
%       % Check to see if either has completed using GOP which
%       % combines the results of 'gotSolution' from each lab
%       done = gop(@any, gotSolution);
%     end
%   end
% access solution in 'state'
A =fopen('dataset.txt','r');
U = fscanf(A,'%f');
dt=0.1;
r_t=0;
max=0.9;
L=10;
l_max=L;
d=0.99;
p=0;
B=3.667;
exceeded =false;
initial=true;
seakeeping =false;
j=1;
n=0;
maxR=0.9;
above_beta=0;
t_end= size(U)/10;
t_vect =0:dt:t_end;

V_vect=zeros(1,length(t_vect));
I_max=1;
I_e_vect =ones(1,length(t_vect));
beta_vect =ones(1,length(t_vect))*8*0.21159;
V_vect(1)=0;

inputlayersize=1;
outputlayersize=1;
hiddenlayersize=7;


%  delete(gcp('nocreate'));
%         parpool(2);
t=cputime;

for i=1:size(U)
    
    if (abs(p)>B)
        exceeded=true; % in the danger zone
        if(abs(p)>L)
            seakeeping=true;
        else
            seakeeping=false;
        end
    else
        exceeded=false;
        if(initial==false)
            initial=true; % next time it is exceeded it is going to be the initial value
        end
    end
    if(seakeeping)
        x=L/B;
        beta=L;
        l_max=abs(p)*x;
    else
        beta=B;
        l_max=L;
    end
    if(~exceeded)% i.e. p<beta
        
        n=d*n;
        if (p>0)
            p= p+U(i)-n;
        else
            p= p+U(i)+n;
        end
        
        [V_vect(i+1),I_e_vect(i),above_beta]= runSNN(abs(p),l_max,V_vect(i),above_beta);
        
    end
    
    if (exceeded)
        if(initial)
            n=initialEstimate(p,U(i),maxR);
            initial=false;
        end
        %         %Run DNN
        
        % To send data from the workers, create a DataQueue object.
        %         Q = parallel.pool.DataQueue;
        %         afterEach(Q,@(data) updateValues(data));
        %         To represent function executions on parallel workers and hold their results, use future objects.
         f(1:3) = parallel.FevalFuture;
        
        f(1)= parfeval(@optimiseParameters, 1, n , beta);
        
        f(2) = parfeval(@runDNN, 1, n, inputlayersize,outputlayersize,hiddenlayersize,false);
        
        f(3)=parfeval(@runDNN, 1, getOutput((abs(p)+U(i)-beta),maxR), inputlayersize,outputlayersize,hiddenlayersize,true);
        results = cell(1,3);
        
        for idx = 1:3
            [completedIdx,value] = fetchNext(f);
            results{completedIdx} = value;
            %fprintf('Got result with index: %d. the value of %d\n', completedIdx, value);
        end
        if(results{1} <= maxR)
            if(results{2}<maxR)
                n=getOutput(results{1},results{2});
            else
                n=results{1};
            end
            
        elseif(results{2}<maxR)
            n=results{2};
        else
            n=maxR;
        end
        
        
      
   

        %        n=d*n+g;
        
        %        T=t+T;
        %%
        
        if (p>0)
            p= p+U(i)-n;
        else
            p= p+U(i)+n;
        end % Calculate the position
        
        
        %%
        [V_vect(i+1),I_e_vect(i),above_beta]= runSNN(abs(p),l_max,V_vect(i),above_beta);
        
        
        if abs(p)<B
            exceeded=false;
        end
        %
    end
    
    
    
end

t=cputime-t

plot (t_vect,V_vect);
hold on;
% plot(t_vect,beta_vect, '-r');
title ('Voltage vs Time');
xlabel('Time in ms')
ylabel('Voltagein mV');


function r_t= initialEstimate(position,env,maxR)
a =4.216;
b =0.5771;
c=0.1804;
m=0.03345;
n= 4.164;
r_p=(a*(position)+b)/(position+c);

r_e=(m*env)+n;
if(r_p<maxR)&&(r_e<maxR)
    r_t=max(r_p,r_e);
elseif(r_p<maxR)
    r_t=r_p;
elseif(r_e<maxR)
    r_t=r_e;
else
    r_t=maxR;
end

end




function n=getOutput(n,max)
if n>max
    n=max;
end
end
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

%indicate spike
if (V_vect>V_th)
    
    V_vect =V_reset;
    
end


end
