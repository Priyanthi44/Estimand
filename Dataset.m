function  Y=Dataset(side,num,items)



switch num
    case 1
        
        Y(1:items,1)=0.9;
%         Y =force.*Y;%constant number
%         fileID = fopen(strcat('z', num2str(1),'.txt'),'w');
%         fprintf(fileID,'Constant\n\n');
%         fclose(fileID);
%         fileID = fopen(strcat( num2str(side),'.txt'),'w');
%         fprintf(fileID,'%2.3f\n',Y);
%         fclose(fileID);
    case 2
        
         Y=rand(items,1);
         m = mean(Y); %uniform distribution
        st= std(Y);
        Y= (Y-m)/st;
%          Y=abs(Y);
%          Y(Y>1) =1;
%          Y =force.*Y;
%          fileID = fopen(strcat('z', num2str(2),'.txt'),'w');
%          fprintf(fileID,'Uniform\n\n');
%          fclose(fileID);
%          fileID = fopen(strcat( num2str(side),'.txt'),'w');
%         fprintf(fileID,'%2.3f\n',Y);
%         fclose(fileID);
    case 3
        
         Y=(randn(items,1));% normal distribution
            m = mean(Y);
        st= std(Y);
        Y= (Y-m)/st;
       
%          Y =force.*Y;
%         fileID = fopen(strcat('z', num2str(3),'.txt'),'w');
%         fprintf(fileID,'Normal\n\n');
%         fclose(fileID);
%         fileID = fopen(strcat( num2str(side),'.txt'),'w');
%         fprintf(fileID,'%2.3f\n',Y);
%         fclose(fileID);
    case 4
        
        Y= zeros(items,1);
%          for i=1:items
%              Y(i,1) =3*sin(i);
%              % Fs=5000; % sampling fre
% 
%          end
         f=200;
n=[0:1/items:1];
Y=2*sin(2*pi*f*n);

%             m = mean(Y); %uniform distribution
%         st= std(Y);
%         Y= (Y-m)/st;
%          Y =force.*Y;
%          fileID = fopen(strcat('z', num2str(4),'.txt'),'w');
%          fprintf(fileID,'Sin wave\n\n');
%          fclose(fileID);
%          fileID = fopen(strcat( num2str(side),'.txt'),'w');
%         fprintf(fileID,'%2.3f\n',Y);
%         fclose(fileID);
    case 5
       
         Y= zeros(items,1);% Linear input
             for i=1:items
                 Y(i,1) =(i);
%                  if Y(i,1)<0
%                      Y(i,1)=abs(Y(i,1));
%                  end
             end
%                 m = mean(Y); %uniform distribution
%         st= std(Y);
%         Y= (Y-m)/st;
%            for i=1:items
%                  
%                  if Y(i,1)<0
%                      Y(i,1)=abs(Y(i,1));
%                  end
%            end
%              Y =force.*Y;
%              fileID = fopen(strcat('z', num2str(5),'.txt'),'w');
%              fprintf(fileID,'Linear\n\n');
%              fclose(fileID);
%              fileID = fopen(strcat( num2str(side),'.txt'),'w');
%         fprintf(fileID,'%2.3f\n',Y);
%         fclose(fileID);
    case 6  % Step Function
       
        start =ceil(items/3);

             Y= zeros(items,1);
             for i=start:items
                 Y(i,1) =0.2;
             end
             
%              Y =force.*Y;
%              fileID = fopen(strcat('z', num2str(side),'.txt'),'w');
%              fprintf(fileID,'Step Function\n\n');
%              fclose(fileID);
%              fileID = fopen(strcat( num2str(6),'.txt'),'w');
%         fprintf(fileID,'%2.3f\n',Y);
%         fclose(fileID);
      case 7 
    
        Y= zeros(items,1); % Exponential curve
         for i=1:items
             Y(i,1) =0.001*exp(i);
         end
% %            m = mean(Y);
% % %         st= std(Y);
%         Y= Y/m;
%          Y =force.*Y;
%          fileID = fopen(strcat('z', num2str(4),'.txt'),'w');
%          fprintf(fileID, 'Exponential Function\n\n');
%          fclose(fileID);
%          
%          fileID = fopen(strcat( num2str(side),'.txt'),'w');
%         fprintf(fileID,'%2.3f\n',Y);
%         fclose(fileID);
    otherwise
            Y= zeros(items,1);
end
           
end
 %persistent Y;
% global l;
% global A;
% global B;
%% Constant input
 % Y(1:items,1)=0.9;
 %% Exponential Input
%  
%  Y= zeros(items,1);
%  for i=1:items
%      Y(i,1) =exp(i)*0.000000001;
%  end
% % %  
 %% Linear
%  Y= zeros(items,1);
%  for i=1:items
%      Y(i,1) =0.000003*(i)+1;
%  end
 
  %% Step
 
%   start =2000;
% 
%  Y= zeros(items,1);
%  for i=start:items
%      Y(i,1) =0.81;
%  end
%  max(Y)
% %% Quadratic
%   Y= zeros(items,1)
%  for i=1:items
%      Y(i,1) =0.000003*(i^2);
%  end

%% Cube
%   Y= zeros(items,1)
%  for i=1:items
%      Y(i,1) =0.0003*(i^3)+1;
%  end
%% Sin
%   Y= zeros(items,1);
%  for i=1:items
%      Y(i,1) =sin(i)+1;
%  end
%% Random Data -uniform
%  Y=rand(items,1);
%% Random Data- normal
%  Y=(randn(items,1));
% X= cat(1,A,B);
 % Y=awgn(Y,10, 'measured');
% max(Y)
% l=mean(Y);
% A=min(Y);
% B=max(Y);
% %  X=X+Y;
 

% Add environmental turbulance
% Fs=5000; % sampling fre
% f=200;
% n=[0:1/Fs:5];
% X=sin(2*pi*f*n);
% plot(X);
% Y = cat(1,Y,readNormalLinebyLine(strcat('normal/' ,int2str(1),'.txt')));
% mean(X)
% X = zeros(2,size(Y,1));
% X=Y;

% function [u1,u2]= readLinebyLine(filename) %, from, to)
%       fid = fopen(filename);
% 
% tline = fgets(fid);
% while ischar(tline)
%     
%    A = strsplit(tline);
%     tline = fgets(fid);
%     %if (str2double(A(1)) >from && str2double(A(1))<to)
%     u1 =str2double(A(27));
%     u2=str2double(A(29));
%   
%    % end
% end
% 
% 
% fclose(fid);
% function X= readNormalLinebyLine(filename) %, from, to)
%       fid = fopen(filename);
% X =[0 0];
% tline = fgets(fid);
% while ischar(tline)
%     
%    A = strsplit(tline);
%     tline = fgets(fid);
%     %if (str2double(A(1)) >from && str2double(A(1))<to)
%     u1 =str2double(A(27));
%     u2=str2double(A(29));
%     X =[ X;horzcat(u1,u2)];
%    % end
% end
% end