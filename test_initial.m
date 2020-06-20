
clear all
close all;
clc;

xm=100;
ym=100;

sink.x=0.5*xm;
sink.y=0.5*ym;

n=100;

for i=1:1:n
    S(i).Initial_energy= 0.5*rand(1,1)+0.5;
    S(i).RE=S(i).Initial_energy;
end

ETX=50e-9;
ERX=50e-9;
%Transmit Amplifier types
Efs=10e-12;              
Emp=0.0013e-12;         
%Data Aggregation Energy
EDA=5e-9;

%maximum number of rounds
rmax = 4000 ;
%Computation of do
d0=sqrt(Efs/Emp);

Tc = 20;

fis1 = readfis('dis_Fuzzyfitness1');           
fis2 = readfis('dis_Fuzzyfitness2');               
fis3 = readfis('Cluster.radius');

%Creation of the random Sensor Network
figure(1);
for i=1:1:10
    for j=1:1:10
        S((i-1)*10+j).xd= ((i-1) + rand(1))*10;
        XR((i-1)*10+j)=S((i-1)*10+j).xd;
        S((i-1)*10+j).yd=((j-1) + rand(1))*10;
        YR((i-1)*10+j)=S((i-1)*10+j).yd;
        S((i-1)*10+j).G=0;
        %initially there are no cluster heads only nodes
        S((i-1)*10+j).type='N';
        % Node identification
        S((i-1)*10+j).id=S((i-1)*10+j);
        S((i-1)*10+j).state = 'Initial_State';
        %Random node distribution
        plot(S((i-1)*10+j).xd,S((i-1)*10+j).yd,'o');
        hold on;
    end
end
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
S(n+1).type  = 'BS';
S(n+1).id = n+1;
S(n+1).RE = 0.5;
max = 0;
min =  10000;
R_comp = 20;
c = 0.4;
bit = 4000;
Rmax = 30;
for i = 1:1:n
    S(i).distobase = sqrt((S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
    if(S(i).distobase > max )
        max = S(i).distobase;
    end
    if(S(i).distobase < min)
        min = S(i).distobase;
    end 
end

for i = 1:1:n
    S(i).R = (1 - c*(max - S(i).distobase)/(max - min))*R_comp;
end

plot(S(n+1).xd,S(n+1).yd,'x');

save('test_init.mat');