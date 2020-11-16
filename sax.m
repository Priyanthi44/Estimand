% startRange = 2;
%  stdc= 1;
%  endRange = 512;
% 
%  table = cell(endRange-startRange,1);
%   for r=startRange:endRange
%     table{r-startRange+1} = norminv((1:r-1)/r,0,stdc);
%   end

%% LOAD DATA and NORMALISE
%DNN data
  A =fopen('nfeedback.txt','r');
    valArray = fscanf(A,'%f');
    zmean=mean(valArray);
    zstd=std(valArray);
    valArray = (valArray-zmean)/zstd;
    
     % Random data
    A =fopen('random.txt','r');
    rdArray = fscanf(A,'%f');
    rdArray = rdArray(size(rdArray)-size(valArray,1)+1:end);
    zmean=mean(rdArray);
    zstd=std(rdArray);
    rdArray = (rdArray-zmean)/zstd;
    
    % Position data
    A =fopen('position.txt','r');
    posArray = fscanf(A,'%f');
    zmean=mean(posArray);
    zstd=std(posArray);
    posArray = (posArray-zmean)/zstd;
    
 %% CREATE PAA  
    
  w =32;
  p=1;
  %q=size(valArray,1);
%   E_x=getEsumx(valArray,q)-getEsumx(valArray,p);
%   E_xx=getEsumxx(valArray,q)-getEsumxx(valArray,p);
  L=size(valArray,1);
%   mu=E_x/L;
%   sigma=sqrt((E_xx-E_x^2/L)/L-1);
mu_val=mean(valArray);
mu_rd =mean(rdArray);
mu_pos=mean(posArray);

sigma_val=std(valArray);
sigma_rd=std(rdArray);
sigma_pos=std(posArray);
  segments=fix(L/w);
  paa_val=zeros(1,segments);
  paa_rd=zeros(1,segments);
  paa_pos=zeros(1,segments);
  start=1;
  p_end=w-1;
  
  for i=1:segments
%       paa(i)=((getEsumx(valArray,p_end)-(getEsumx(valArray,start))/(L/w))-mu)/sigma;
paa_val(i)=(getSum(start,p_end,valArray)/w-mu_val)/sigma_val;
paa_rd(i)=(getSum(start,p_end,rdArray)/w-mu_rd)/sigma_rd;
paa_pos(i)=(getSum(start,p_end,posArray)/w-mu_pos)/sigma_pos;
      start=p_end+1;
      p_end=start+w-1;
  end
%   max(paa_val)
%   min(paa_val)
%     max(paa_rd)
%   min(paa_rd)
%     max(paa_pos)
%     min(paa_pos)
%   figure(1);
%        [counts, bins] = hist(paa);
%    plot(bins, counts);
%    t_vect=1:1:size(paa,2);
%    figure(2);
%  plot(t_vect,paa);

%% CREATE SAX
isax_val= convertPaatoSax(paa_val);
isax_rd= convertPaatoSax(paa_rd);
isax_pos= convertPaatoSax(paa_pos);

isax_val=char(isax_val);
isax_rd=char(isax_rd);
isax_pos=char(isax_pos);

wordlen=8;

valLinkedlist= strings(fix(size(isax_val,2)/wordlen)-1,1);
rdLinkedlist= strings(fix(size(isax_rd,2)/wordlen)-1,1);
posLinkedlist= strings(fix(size(isax_pos,2)/wordlen)-1,1);

for i=1:fix(size(isax_val,2)/wordlen)-2
    valLinkedlist(i)=string(isax_val(i:i+wordlen-1));
    rdLinkedlist(i)=string(isax_rd(i:i+wordlen-1));
    posLinkedlist(i)=string(isax_pos(i:i+wordlen-1));
end

valLinkedlist=categorical(valLinkedlist);
posLinkedlist=categorical(posLinkedlist);
rdLinkedlist=categorical(rdLinkedlist);

%% IDENTIFY MOTIFS OF DNN
figure
h1 = histogram(valLinkedlist);
 xlabel("Words")
 ylabel("Frequency")
title("DNN Distribution")
valCounts = h1.BinCounts;
valNames = h1.Categories;

idxCounts = valCounts > 3;
Mofits = valNames(idxCounts);
idxMotifs = ismember(valLinkedlist,Motifs);
data(idxInfrequent,:) = [];


idxMaxCounts = valCounts == max(valCounts);
frequentClasses = valNames(idxMaxCounts);
idxfrequent = ismember(valLinkedlist,frequentClasses);

count=0;
motifpos=zeros(max(valCounts),1);

for i=1:size(idxfrequent)
    if idxfrequent(i)==1
        count=count+1;
        motifpos(count)=i;     
    end
    if(count==max(valCounts))
        break;
    end
end
%% EXPAND THE MOTIF
%move forward the motif until a repeating pattern is not found    
 wordlenbefore = wordlen-compareElements(valLinkedlist(motifpos-1),0)
 wordlenafter = compareElements(valLinkedlist(motifpos+1),1)
 if (wordlenbefore==wordlen)
     wordlenbefore=0;
 end
 

%% COMPARE SUBSEQUENCES

valBefore=char(valLinkedlist(motifpos-1));
valAfter=char(valLinkedlist(motifpos+1));
valList =strcat(valBefore(wordlen-wordlenbefore:wordlen),char(valLinkedlist(motifpos)),valAfter(1:wordlenafter));

rdBefore=char(rdLinkedlist(motifpos-1));
rdAfter=char(rdLinkedlist(motifpos+1));
rdList =strcat(rdBefore(wordlen-wordlenbefore:wordlen),char(rdLinkedlist(motifpos)),rdAfter(1:wordlenafter));

posBefore=char(posLinkedlist(motifpos-1));
posAfter=char(posLinkedlist(motifpos+1));
posList=strcat(posBefore(wordlen-wordlenbefore:wordlen),char(posLinkedlist(motifpos)),posAfter(1:wordlenafter)) ;

%% TODO ompare rdList items with each other  and posList items with each other and identify motifs if any  
rdMap=bitMap(rdList);
posMap=bitMap(posList);
valMap=bitMap(valList);
image(rdMap,'CDataMapping','scaled');
colorbar
image(posMap,'CDataMapping','scaled');
colorbar
image(valMap,'CDataMapping','scaled');
colorbar
% pos =1:1:10;
% figure(2);
% plot(pos,motifpos);
function k= compareElements(List,order)
match=true;
len=strlength(char(List(1)));

if(order)
    k=1;
    
while(k<=len)
for j=1:size(List)-1
    
    el1=List(1);
    el2=List(j+1);
    
    
         el1=char(el1);
         el2=char(el2);
        
             if(el1(k)==el2(k))
                 match=true;
             else
                 match=false;
                 
             end
            
        
 
 
    if(~match)
        break;
    end
    
end
  if(~match)
      k=k-1;
      break;
        
   end
         k=k+1;
     
   
end


else
     k=len;
  
 while(k>=1)
    
for j=1:size(List)-1
    
    el1=List(1);
    el2=List(j+1);

    
    
         el1=char(el1);
         el2=char(el2);
        
             if(el1(k)==el2(k))
                 match=true;
             else
                 match=false;
             end
            
        
 
 
    if(~match)
        break;
    end
     
end
    if(~match)
        break;
    end
    k=k-1;
    
  
 end 
end

end

function isax= convertPaatoSax(paa)
isax=zeros(1,size(paa,2));
r=6;
range=norminv((1:r-1)/r,0,1);
for i=1:size(paa,2)
        if paa(i)>range(5)
        isax(i)='f';
        elseif paa(i)>range(4)
        isax(i)='e';
        elseif paa(i)>range(3)
        isax(i)='d';
        elseif paa(i)>range(2)
        isax(i)='c';
        elseif paa(i)>range(1)
        isax(i)='b';
        else 
        isax(i)='a';
        end
        
end
end

function bmap= bitMap(L)
bmap=zeros(6,6);
s=size(L);
List=L';
len=s(1)*s(2);
nexteven=0;

for i=1:len
    
    if(mod(i,23))
        if(i>23)
            if(nexteven)&&(~mod(i,2))
               el=strcat(List(i),List(i+1));
               x=getIndex(el);
               bmap(x(1),x(2))=bmap(x(1),x(2))+1; 
            elseif(mod(i,2))   
                el=strcat(List(i),List(i+1));
                x=getIndex(el);
                bmap(x(1),x(2))=bmap(x(1),x(2))+1;
            end  
        
        elseif(mod(i,2))   
        el=strcat(List(i),List(i+1));
        x=getIndex(el);
        bmap(x(1),x(2))=bmap(x(1),x(2))+1;
        end
      
    else
        if(mod(i/23,2))
            nexteven=1;
            
        else
            nexteven=0;
            
        end
        
    end
        
end


end

function indx=getIndex(el)
indx= [0,0];
switch el
    case 'aa'
        indx=[1,1];
    case 'ab'
        indx=[1,2];
    case 'ac'
        indx= [1,3];
    case 'ad'
        indx=[1,4];
    case 'ae'
        indx=[1,5];
    case 'af'
        indx= [1,6];
     case 'ba'
        indx=[2,1];
    case 'bb'
        indx=[2,2];
    case 'bc'
        indx= [2,3];
    case 'bd'
        indx=[2,4];
    case 'be'
        indx=[2,5];
    case 'bf'
        indx= [2,6];
    case 'ca'
        indx=[3,1];
    case 'cb'
        indx=[3,2];
    case 'cc'
        indx= [3,3];
    case 'cd'
        indx=[3,4];
    case 'ce'
        indx=[3,5];
    case 'cf'
        indx= [3,6];  
          case 'da'
        indx=[4,1];
    case 'db'
        indx=[4,2];
    case 'dc'
        indx= [4,3];
    case 'dd'
        indx=[4,4];
    case 'de'
        indx=[4,5];
    case 'df'
        indx= [4,6];
    case 'ea'
        indx=[5,1];
    case 'eb'
        indx=[5,2];
    case 'ec'
        indx= [5,3];
    case 'ed'
        indx=[5,4];
    case 'ee'
        indx=[5,5];
    case 'ef'
        indx= [5,6];
      case 'fa'
        indx=[6,1];
    case 'fb'
        indx=[6,2];
    case 'fc'
        indx= [6,3];
    case 'fd'
        indx=[6,4];
    case 'fe'
        indx=[6,5];
    case 'ff'
        indx= [6,6];
    
end
end
 function sum=getSum(istart,iend, array)
 sum=0;
 for i=istart:iend
     sum=sum+ array(i);
 end
 end
  function Esum_x =getEsumx(valArray,iend)
  Esum_x=0;
  for i=1:iend
      Esum_x=Esum_x+valArray(i);
  end
  end 
  function Esum_xx =getEsumxx(valArray,iend)
  Esum_xx=0;
  for i=1:iend
      Esum_xx=Esum_xx+valArray(i)^2;
  end
  end
  
 
  