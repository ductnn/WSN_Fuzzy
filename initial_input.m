xm=100;
ym=100;
 
%x and y Coordinates of the sink
sink.x=50;
sink.y=50;
 
%Number of Nodes in the Field
n=100;
 
%Probability of a node to become cluster head
%p=0.1;
Tc= 20;
 
Rmax = 30;
 
%Eelc=Etx=Erx
ETX=50* 10^-9;
ERX=50* 10^-9;
%Transmit Amplifier types
Efs=10* 10^-12;              
Emp=0.0013* 10^-12;            
%Data Aggregation Energy
EDA=5* 10^-9;
 
%maximum number of round 
rmax = 1000;
 
%Computation of do
d0 = sqrt(Efs/Emp);
 
bit = 4000;
fis1 = readfis('dis_Fuzzyfitness1');           
% fis2 = readfis('dis_Fuzzyfitness2');               
% fis3 = readfis('Cluster.radius');

%Creation of the random Sensor Network
figure(1);

for i=1:1:n
    %Initial Energy
    S(i).Initial_energy = 0.5;
    S(i).RE = S(i).Initial_energy;
    S(i).xd=rand(1,1)*xm;
    S(i).yd=rand(1,1)*ym;
    %initially there are no cluster heads only nodes
    S(i).type = 'N';
    S(i).id=i;
    S(i).state = 'Initial_State';
    S(i).distoBS = norm([S(i).xd-50 S(i).yd-50]);
    plot(S(i).xd,S(i).yd,'o','MarkerFaceColor','b');
    hold on;
end
S(n+1).xd = sink.x;
S(n+1).yd = sink.y;
plot(sink.x, sink.y, 's', 'MarkerFaceColor', 'red');
text(sink.x, sink.y,'  BS','Color','b','FontWeight','b');
hold on;
grid on;
save('fibo_test_1.mat');