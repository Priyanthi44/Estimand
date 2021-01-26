% In this function, we attempt to optimise the parameters of the equations
%in here we apply Nelder-Mead algorithm to optimise parameters
%need persistent parameters
%if the beta<position, the initial N is calculated using the equations
function r_t=optimiseParameters(r_t,beta)
persistent a
if isempty(a)
    a=0;
    
end
persistent b

if isempty(b)
    b=0;
    
end
persistent c
if isempty(c)
    c=0;
    
end
persistent m
if isempty(m)
    m=0;
    
end
persistent n
if isempty(n)
    n=0;
    
end
persistent firsttime
if isempty(firsttime)
    firsttime=true;
end

if(firsttime)
    firsttime=false;
%     a =a+4.216;
%     b =b+0.5771;
%     c=c+0.1804;
%     m=m+0.03345;
%     n= n+4.164;
end
%
%  FUN can be a parameterized function. Use an anonymous function to
%    %     capture the problem-dependent parameters:
%    %        f = @(x,c) x(1).^2+c.*x(2).^2;  % The parameterized function.
%    %        c = 1.5;                        % The parameter.
%    %        X = fminsearch(@(x) f(x,c),[0.3;1])
%    %
%   %   FMINSEARCH uses the Nelder-Mead simplex (direct search) method.

%equation
%delta_r= (c *r_t-b)/(a-r_t) +(r_t-n)/m -(beta +r_t);

f=@(x) (x(3) *r_t  -x(2))/ (x(1)-r_t) + (r_t- x(5))/x(4) -(beta+r_t)  ;
X=fminsearch(f,[a;b;c;m;n]);
if(X(4)<X(5))
    f=@(x) (x(3) *r_t  -x(2))/ (x(1)-r_t) + (r_t- x(49)) -(beta+r_t)  ;
end
if~(X(4)<=0.01)||~(X(1)>=r_t)
    a= X(1);
    b= X(2);
    c= X(3);
    
    m= X(4);
    n= X(5);
end
r_t=r_t+ (c *r_t-b)/(a-r_t) +(r_t-n)/m -(beta +r_t);



end


