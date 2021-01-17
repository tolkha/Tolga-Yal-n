function[xk,i]=Goldensec(F,xl,xh,xend)
syms x1 
tol=xend/(xh-xl);
ni=-2.078*log(tol); %Number of iterations
T=0.38197;

a1=xl+T*(xh-xl);
f1=vpa(subs(F,x1,a1));
a2=xh-T*(xh-xl);
f2=vpa(subs(F,x1,a2));
i=1;
while 1
    if i<ni
        if f1>f2
            xl=a1;
            a1=a2;
            f1=f2;
            a2=xh-T*(xh-xl);
            f2=vpa(subs(F,x1,a2));
        elseif f1<f2
                xh=a2;
                a2=a1;
                f2=f1;
                a1=xl+T*(xh-xl);
                f1=vpa(subs(F,x1,a1));
        end
    else
        xk=(a1+a2)/2;
          fprintf('<GOLDEN SECTÝON>The Number of Iteration is=> %.4f\n',i);
          fprintf('<GOLDEN SECTÝON>The Root is=%.4f\n',xk);
          fprintf('\n');
        break
    end
    i=i+1;
end
end