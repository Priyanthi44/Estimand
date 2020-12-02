A =fopen('input.txt','r');
U = fscanf(A,'%f');

num1 =0;
max(U)
for i=1:size(U)
    total = 2020-U(i);%first number U(i)
    
    for j=1:size(U)
        
        
        num1 =total-U(j);
        
        %        if (any(U(i:end)) == num)
        %            sum=num+U(i)
        %             answer = num *U(i)
        %             break;
        %        end
        if(num1>0)
            if(num1>min(U))
                if(ismember(num1,U))
                    sum=num1+U(j)+U(i);
                    if(sum==2020)
                        num1
                        answer =num1* U(j) *U(i)
                        break;
                    end
                end
            end
        end
        
        
    end
    
end