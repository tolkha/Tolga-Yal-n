function [yi]=Test_ANN(X,test,s)
    syms x
    [r,n]=size(test);
    
    %wg,bh,wc,bc Calculation
    for i=1:r
            wg(:,i)=X((s*i)-(s-1):(s*i));
    end
    bh=X(s*r+1:s*r+s,1);
    wc(1,:)=X(s*(r+1)+1:s*(r+2),1);
    bc=X(s*(r+2)+1);
    %wg,bh,wc,bc Calculation
     
    h=(exp(x)-exp(-x))/(exp(x)+exp(-x));
        
    %yi Calculation
    for i=1:n
            yi(i,1)=wc*subs(h,x,wg*test(:,i)+bh)+bc;
    end
    %yi Calculation
end