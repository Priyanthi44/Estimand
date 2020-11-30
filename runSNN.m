% I know the number of point means x axis is known that is t_vect
% I want the graph to be plotted once the simulation finishes I want to plot
% V_vect and t_vect

function runSNN(i,iend)
dt=0.1;
R_m =10;
tau =100;
Th=8;
I_max=1;
V_th=Th;
V_reset = 0;
% V_vect and I_vect remembers the value
persistent  V_vect;
persistent I_vect;
t_vect =0:dt:t_end;
if(isEmpty(I_vect))
    I_vect =ones(1,length(t_vect));
end
if(isEmpty(V_vect))
    V_vect =zeros(1,length(t_vect));
end
err=abs((p/l)*100); % Percentage of error
I_vect(i)=I_max*(1-err/100); % Adjust input current according to the error
if(err>100)
    I_vect(i)=0;
end
V_inf = I_vect(i)* R_m;
V_vect(i+1) =V_inf + (V_vect(i)-V_inf) *exp(-dt/tau); %Calculate voltage
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

%indicate spike
if (V_vect(i+1)>V_th)
    
    V_vect(i+1) =V_reset;
    
end

%if i==intems, then plot the graph
if(i==iend)
    figure;
    t_vect =0:dt:iend/10;
    beta_vect =ones(1,length(t_vect))*8*0.21159; 
    plot (t_vect,V_vect);
    hold on;
    plot(t_vect,beta_vect, '-r');
    title ('Voltage vs Time');
    xlabel('Time in ms')
    ylabel('Voltagein mV');
    
end


end