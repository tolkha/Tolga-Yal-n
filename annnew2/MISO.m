function [X,o,yi_value,error_E]=MISO(s,r,t,n,MSE,o_max,yr)
        P=rand(s*(r+2)+1,1); % Starting Points for L_M
        X(:,1)=P; 
        syms x [1 length(P)];
        x=x.';
        
        h=(exp(x(1))-exp(-x(1)))/(exp(x(1))+exp(-x(1))); %Activation function

        n_max=100; %L_M için Maksimum Ýteraston Sayýaý
        mu=1; %L_M için mü (Between mü_min and mü_max) 
        mu_scale=10; %L_M için mü scale (Between 0 and 1)
        mu_min=0.0001; %L_M için mü minimum
        mu_max=10000; %L_M için mü maximum
        E1=10^-4; %L_M için Koþul 2
        E2=10^-4; %L_M için Koþul 3
        E3=10^-4; %L_M için Koþul 4
        b=0;

        %wg,bh,wc,bc Calculation
        for i=1:r
            wg(:,i)=x(s*i-(s-1):s*i);
        end
        bh=x(s*r+1:s*r+s,1);
        wc(1,:)=x(s*(r+1)+1:s*(r+2),1);
        bc=x(s*(r+2)+1);
        %wg,bh,wc,bc Calculation

        %yi Calculation
        for i=1:n
            yi(i,1)=wc*subs(h,x1,wg*t(:,i)+bh)+bc;
        end
        %yi Calculation

        %Jacobian Calculation
        for i=1:length(t)
            for j=1:s*r
                k=mod(j-1,s)+1;
                m=fix((j-1)/s)+1;
                h2_value=vpa(subs(h,x1,wg(k,:)*t(:,i)+bh(k,1)));
                J(i,j)=-(wc(k)*t(m,i)*(1-h2_value^2));
            end
        end
        for i=1:length(t)
            for j=s*r+1:s*r+s
                h2_value=vpa(subs(h,x1,wg(j-s*r,:)*t(:,i)+bh(j-s*r,1)));
                J(i,j)=-(wc(j-s*r)*(1-h2_value^2));
            end
        end
        for i=1:length(t)
            for j=s*r+s+1:s*r+2*s
                J(i,j)=-(vpa(subs(h,x1,wg(j-(r+1)*s,:)*t(:,i)+bh(j-(r+1)*s,1))));
            end
        end
        for i=1:length(t)
            for j=s*r+2*s+1
                J(i,j)=-1;  
            end
        end
        %Jacobian Calculation

        E=vpa((yr-yi),4);
        F=vpa(E.'*E,4);

        o=1;
        while 1
            [P,i]=L_M(F,P,E,J,mu,mu_scale,mu_min,mu_max,n_max,E1,E2,E3);
            fprintf('Number of iteration of MISO is %d\n',o);
            fprintf('Number of neuron is %d\n',s);
            P=P(:,i+1);
            X(:,o+1)=P;
            error_E(o,1)=norm(subs(E,x(:,1),X(:,o+1)));
            delta_X=norm(X(:,o+1)-X(:,o));
            if delta_X==0
                fprintf('L_M is in a loop please change your parameters\n');
                b=b+1;
            end
            if error_E(o,1)<MSE
                fprintf('Target error has been reached.\n');
                yi_value=vpa(subs(yi,x(:,1),X(:,o+1)));
                disp('    yi Values');
                disp(yi_value);
                a=0;
                break;
            elseif o>o_max
                fprintf('Maximum number of iterations is exceeded.\n');
                yi_value=double(subs(yi,x(:,1),X(:,o+1)));
                disp('    yi Values');
                disp(yi_value);
                break;
            end
            o=o+1;
        end
end