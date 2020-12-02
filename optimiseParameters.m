% In this function, we attempt to optimise the parameters of the equations
%in here we apply Nelder-Mead algorithm to optimise parameters
%need persistent parameters
%if the beta<position, the initial N is calculated using the equations
function r_t=optimiseParameters(initial,position,env)
persistent a;
persistent b;
persistent c;
persistent m;
persistent n;
if(isEmpty(a))
    a =0;
end
if(isEmpty(b))
    b =0;
end
if(isEmpty(c))
    c =0;
end
if(isEmpty(m))
    m=0;
end
if(isEmpty(n))
    n =0;
end
r_t=0;
if(initial)
    r_t= initialEstimate(position, env,a,b,c,m,n);
else
    
end

end

function r_t= initialEstimate(a,b,c,m,n)
r_t=0;
r_t=
end