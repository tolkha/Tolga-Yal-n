clc;
clear;
syms x

s=10; %Nöron Sayýsý
o_max=100; %Maximum iteration number for MISO
MSE=10^-1;

data=xlsread('databaseexport13122018.xlsx'); %Taking datas from Excel file
t= data(1:45,3:4); % MISO inputs
t=t.';
[r,n]=size(t);
yr=data(1:n,7); %Real output

[X,o,yi_value,error_E]=MISO(s,r,t,n,MSE,o_max,yr);
%plot(yr);
%hold on;
%plot(yi_value,'r*');

test=data(45:end,3:4);
test=test.';
[yi]=Test_ANN(X(:,o+1),test,s);
datasonuc=xlsread('databaseexport13122018.xlsx'); %Taking datas from Excel file
plot(yi,'r*');
 
t1= data(45:end,7);% MISO inputs 
hold on 
plot(t1)
for i=1:size(t1)
RMSE(i) = sqrt(mean(t1(i)-yi(i)).^2)
end
