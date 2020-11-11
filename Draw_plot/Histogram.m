load('fibo_test_1.mat');
figure(1);
S(n+1).xd = 50;
S(n+1).yd = 50;
a = [S.xd];
b= [S.yd];
h1=plot(a(1:1:100),b(1:1:100),'bo','MarkerFaceColor','b','DisplayName','Sensor');
hold on
h2=plot(a(101),b(101),'r^','MarkerFaceColor','r','DisplayName','BS');
grid on;
h = [h1(1);h2];
legend(h,'Sensor','BS');
set(gca,'GridLineStyle','--');