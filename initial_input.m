xm=100;
ym=100;
 
%x and y Coordinates of the sink
sink.x=0.5*xm;
sink.y=0.5*ym;

 
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
fis2 = readfis('dis_Fuzzyfitness2');               
fis3 = readfis('Cluster.radius');

%Creation of the random Sensor Network
figure(1);
count = 0;

for i=1:1:10
    for j=1:1:10
        S((i-1)*10+j).xd= ((i-1) + rand(1))*10;
%         XR((i-1)*10+j)=S((i-1)*10+j).xd;
        S((i-1)*10+j).yd=((j-1) + rand(1))*10;
%         YR((i-1)*10+j)=S((i-1)*10+j).yd;
%         S((i-1)*10+j).G=0;
        %initially there are no cluster heads only nodes
        S((i-1)*10+j).type='N';
        % Node identification
        S((i-1)*10+j).id=(i-1)*10+j;
        S((i-1)*10+j).state = 'Initial_State';
        S((i-1)*10+j).Initial_energy = 0.5;
        S((i-1)*10+j).RE = 0.5;
        %Random node distribution
        plot(S((i-1)*10+j).xd,S((i-1)*10+j).yd,'o');
        hold on;
    end
end

% for j=1:1:30
%     %Initial Energy
%     S(j).Initial_energy = 0.5;
%     S(j).RE = S(j).Initial_energy;
% 
%     % checkbox
%     S(j).xd=rand(1,1)*50 + 20;
%     XR(j)=S(j).xd;
%     S(j).yd=rand(1,1)*50 + 20;
%     YR(j)=S(j).yd;
%     %S(i).G=0;
%     
%     dist_to_BS = sqrt( (S(j).xd-sink.x)^2 + (S(j).yd-sink.y)^2 );
%     if dist_to_BS <= Rmax
%         count = count + 1;
%     end
% 
%     %initially there are no cluster heads only nodes
%     S(j).type = 'N';
%     % Node identification
%     S(j).id=j;
%     S(j).state = 'Initial_State'
%     plot(S(j).xd,S(j).yd,'o');
%     hold on;
% end
% 
% for i=1:1:70
%     %Initial Energy
%     S(i).Initial_energy = 0.5;
%     S(i).RE = S(i).Initial_energy;
%     
%     
%     % checkbox
%     S(i).xd=rand(1,1)*xm;
%     XR(i)=S(i).xd;
%     S(i).yd=rand(1,1)*ym;
%     YR(i)=S(i).yd;
%     %S(i).G=0;
%  
%     %initially there are no cluster heads only nodes
%     S(i).type = 'N';
%     % Node identification
%     S(i).id=i;
%     S(i).state = 'Initial_State'
%  
%     %Random node distribution
%     plot(S(i).xd,S(i).yd,'o');
%     hold on;
% end
 
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'x');

save('test.mat');







 
