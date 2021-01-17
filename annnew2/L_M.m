function [P,i]=L_M(F,P,E,J,mu,mu_scale,mu_min,mu_max,n_max,E1,E2,E3)
    [n,~]=size(P); %Fonksiyondaki de�i�ken say�s�, ba�lang�� noktas� olarak verilen fonksiyonun eleman say�s� ile tespit edildi
    syms x [1 n];%Fonksiyondakide�i�ken say�s�na g�re de�i�kenler olu�turuldu
    [~,c]=size(J);
    grad_F=2.*J.'*E; %Gradient of function
    i=1;
    while 1
        %Step2
        f_value=vpa(subs(F,x(1,:).',P(:,i))); %F(xk) parametrik olarak hesapland�
        J_value=vpa(subs(J,x(1,:).',P(:,i))); %J(xk) parametrik olarak hesapland�
        E_value=vpa(subs(E,x(1,:).',P(:,i))); %E(xk) parametrik olarak hesapland�
        %Step3
        while 1
            zk(:,i)=-(inv(J_value.'*J_value+mu.*eye(c)))*J_value.'*E_value; %Her iterasyonda hesaplanan zk, zk matrisinin s�tunlar�na kaydedilir
            f_xkzk=vpa(subs(F,x(1,:).',P(:,i)+zk(:,i)));
            if f_xkzk<f_value
                PK(:,i)=zk(:,i); %�lerleme y�n� aday� zk kesin ilerleme y�n� pk olarak atand�
                f_sk=vpa(subs(F,x(1,:).',P(:,i)+x(1,1).*PK(:,i)));%f(xk+sk*pk) hesapland�
                [sk,~]=Goldensec(f_sk,-20,20,10^-12); %f(xk+sk*pk) y� minimum yapan sk de�eri hesapland�
                P(:,i+1)=P(:,i)+sk.*PK(:,i); %X(k+1) hesapland�
                mu=mu/mu_scale; %mu(k+1) hesapland�
                break; %Step5'e git
            %Step4
            else
                mu=mu*mu_scale;
                if mu<mu_max && mu>mu_min
                    continue; %Ad�m3'e git
                else
                    break;
                end
            end
        end
        %Step5
        %Bitirme �artlar� hesaplan�r ve kontrol edilir
        if mu<=mu_min || mu>=mu_max
            clc;
            fprintf('Condition 5 established\n');
            P(:,i+1)=P(:,i);
            break;
        end
        
        f1=subs(F,x(1,:).',P(:,i));
        f2=subs(F,x(1,:).',P(:,i+1));
        normgrad=norm(subs(grad_F,x(1,:).',P(:,i+1)));
        error=abs(f2-f1);
        delta_x=abs(norm(P(:,i+1)-P(:,i)));
        if mod(i,4)==0
            clc;
        end
        if i>n_max
            clc;
            fprintf('Condition 1 established\n');
            break
        elseif error<E1
            clc;
            fprintf('Condition 2 established\n');
            break
        elseif delta_x<E2
            clc;
            fprintf('Condition 3 established\n');
            break
        elseif normgrad<E3
            clc;
            fprintf('Condition 4 established\n');
            break
        end
        i=i+1;
    end
end