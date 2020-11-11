clear all;
close all;
clc;
figure(1);
load('mazumdar_case1.mat');
y = DEAD(1:60:length(DEAD));
x = (1:60:length(DEAD));
h1 = plot(x, 100- y, 'bo-','MarkerFaceColor','blue');
set(h1,'LineWidth',2);
hold on;

load('mazumdar_case2.mat');
y = DEAD(1:60:length(DEAD));
x = (1:60:length(DEAD));
y(length(y)+1) = 100;
x(length(x)+1) = x(length(x)) + 1;
h2 = plot(x,100-y,'r^-','MarkerFaceColor','red');
set(h2,'LineWidth',2);
ylabel('Number of alive sensor node');
xlabel('Rounds');
grid on;
hold off;
legend('DFCR','Fibonacci');

set(gca,'GridLineStyle','--');
ylabel('Number of alive sensor node','FontWeight','bold','FontAngle','italic');
xlabel('Rounds','FontWeight','bold','FontAngle','italic');
set(gca,'FontWeight','bold');
