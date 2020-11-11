close all;
clear all;
clc;

load('mazumdar_case2.mat');
DFCR_2 = DEAD;
load('mazumdar_case1.mat');
DFCR_1 = DEAD;
load('Fibo_case1.mat');
Fibonacci_1 = DEAD;
load ('Fibo_case2.mat');
Fibonacci_2 = DEAD;

y1 = [find(DFCR_1,1) find(Fibonacci_1,1); find(DFCR_2,1) find(Fibonacci_2,1)];
figure(1);
hb = bar(y1);
set(gca,'xticklabel',{'scenario1','scenario2'},'FontWeight','bold');
ylabel('Rounds','FontWeight','bold');
legend('DFCR','Fibonacci');


figure(2);
a1=find(DFCR_1==50);
a3=find(Fibonacci_1==50); 
a4=find(DFCR_2==50); 
a6=find(Fibonacci_2==50);
y2 = [a1(1) a3(1);a4(1) a6(1)];
hb = bar(y2);
set(gca,'xticklabel',{'scenario1','scenario2'},'FontWeight','bold');
ylabel('Rounds','FontWeight','bold');
legend('DFCR','Fibonacci');
figure(3);
y2 = [length(DFCR_1) length(Fibonacci_1); 1118 length(Fibonacci_2)];
hb = bar(y2);
set(gca,'xticklabel',{'scenario1','scenario2'},'FontWeight','bold');
ylabel('Rounds','FontWeight','bold');
legend('DFCR','Fibonacci');
