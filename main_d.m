clear;
%%%%%%%%%%%%%%%%%%%%    INITIAL PRAMETERS            %%%%%%%%%%%%%%%%%%%
%Field Dimensions - x and y maximun (in meters)
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
 
 
Rmax = 20;
 
%Initial Energy
for i=1:1:n
    S(i).Initial_energy=0.5;
    S(i).RE=S(i).Initial_energy;
end
 
 
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
 %Creation of the random Sensor Network
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    %S(i).G=0;
 
    %initially there are no cluster heads only nodes
    S(i).type = 'N';
    % Node identification
    S(i).id=i;
    S(i).state = 'Initial_State'
 
    %Random node distribution
    plot(S(i).xd,S(i).yd,'o');
    hold on;
end
 
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'x');
 
for r= 1:1:rmax
% Compute Neigbor desity & neighbor cost
for i= 1:1:n
    number_neighbor_i = 0;
    sigma_neigh_cost = 0;
    for j= 1:1:n
        disJtoI = sqrt((S(i).xd-S(j).xd)^2 + (S(i).yd-S(j).yd)^2);
        if (disJtoI <= Rmax && disJtoI > 0)
            number_neighbor_i = number_neighbor_i + 1;
            sigma_neigh_cost = sigma_neigh_cost + (S(i).xd-S(j).xd)^2 + (S(i).yd-S(j).yd)^2;
        end
    end
    S(i).neigh_des = number_neighbor_i/n; 
    S(i).neigh_cost = (sqrt(sigma_neigh_cost/number_neighbor_i))/Rmax;
end
 
%  CH %
for i= 1:1:n
   if (S(i).RE >0 )
       S(i).distoBS = sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2);
       
       % compute Td 
       fis1 = readfis('dis_Fuzzyfitness1');
       Energy_level = S(i).RE/S(i).Initial_energy;
       S(i).Fuzzy_fitness1 = evalfis([Energy_level S(i).distoBS], fis1);
               
       fis2 = readfis('dis_Fuzzyfitness2');
       S(i).Fuzzy_fitness2 = evalfis([S(i).neigh_des S(i).neigh_cost], fis2);
              
       % random alpha from [0.9 1]
       alpha = rand(1,1) / 10 + 0.9;
       
       S(i).Td = alpha * (1 - S(i).Fuzzy_fitness1) * Tc;
       %%S(i).rad = evalfis([S(i).RE S(i).distoBS S(i).Fuzzy_fitness1], fis1);
       
       %Initial candidates
       S(i).candidate = [];
   end
end
 
% Start Bau CH
% Bau CH
number_normal_node = 100;
cluster = 0;
while number_normal_node > 0
 for i= 1:1:n
   CH_selection = 0;
   min_Td = Tc;
 
   for j= 1:1:n
       % Chon thang co Td nho nhat ma chua phai la CH va khong phai la
       % worker cua thang khac (nam trong ban kinh cua 1 CH nao do)
     if (S(j).RE > 0 && S(j).Td < min_Td && isequal(S(j).type,'N') && length(S(j).candidate) == 0)
        min_Td = S(j).Td;
        CH_selection = j; 
     end
   end
 
   if (i == CH_selection)
       S(i).type = 'CH';
       fis3 = readfis('Cluster.radius');
       rad = evalfis([S(i).Fuzzy_fitness1 S(i).Fuzzy_fitness2], fis3);
       S(i).rad = rad;
       S(i).number_worker = 0;
       cluster = cluster +1;
       C(cluster).xd=S(i).xd;
       C(cluster).yd=S(i).yd;
       plot(S(i).xd,S(i).yd,'k*');
              % compute node j received from i
       for t= 1:1:n
         if ((isequal(S(t).type,'N') || isequal(S(t).type,'W'))&& (S(t).RE >0))
            disJToI = sqrt( (S(i).xd-S(t).xd)^2 + (S(i).yd-S(t).yd)^2 );
            if (disJToI <= S(i).rad)
                k = length(S(t).candidate) + 1;
                S(t).type = 'W';
                S(t).candidate(k) = i;
            end
         end  
       end
   end
 end
 
 number_normal_node = 0;
 for i= 1:1:n
    if (isequal(S(i).type,'N') && (S(i).RE >0))
        number_normal_node = number_normal_node + 1;
    end
 end
end
 
% Chon CH cho worker
for j= 1:1:n
    if (isequal(S(j).type,'W') && (S(j).RE >0) && length(S(j).candidate) > 0)
        candidate = S(j).candidate;
        CH_i = candidate(1);
        dist_Sj_CH_i = sqrt((S(j).xd-S(CH_i).xd)^2 + (S(j).yd-S(CH_i).yd)^2 );
        CH_cost_1 = dist_Sj_CH_i * S(CH_i).distoBS/S(CH_i).RE;
        CH = CH_i; 
        
        for i= 2:1:length(candidate)
            CH_i = candidate(i);
            dist_Sj_CH_i = sqrt((S(j).xd-S(CH_i).xd)^2 + (S(j).yd-S(CH_i).yd)^2 );
            CH_cost = dist_Sj_CH_i * S(CH_i).distoBS/S(CH_i).RE;
            
            if (CH_cost < CH_cost_1) 
                CH = CH_i;
                CH_cost_1 = CH_cost;
            end
        end
        
        S(j).CH = CH;
        S(CH).number_worker = S(CH).number_worker + 1;
    end
end
 
 
for i=1:1:n
    if (isequal(S(i).type,'CH') && S(i).RE > 0)
        S(i).L = S(i).distoBS/Rmax;
        S(i).PN = [];
        S(i).parent = n+1;
    end
end
% --- End Bau CH
 
% Routing
for i=1:1:n
    if (isequal(S(i).type,'CH') && S(i).RE > 0)
        S(i).L = S(i).distoBS/Rmax;
    end
end
 
for i=1:1:n
    if (isequal(S(i).type,'CH') && S(i).RE > 0)
        range = k*Rmax;
        for j=1:1:n
            disJToI = sqrt((S(j).xd-S(i).xd)^2 + (S(j).yd-S(i).yd)^2 );
            if (disJToI <= range && isequal(S(j).type,'CH') && S(j).RE > 0 && S(j).L < S(i).L)
                k = length(S(i).PN) + 1;
                S(i).PN(k) = j;
            end 
        end
    end
end
 
for i=1:1:n
    if (isequal(S(i).type,'CH') && S(i).RE > 0 && length(S(i).PN) > 0)
        PN = S(i).PN;
        PN_j = PN(1);
        parent = PN_j;
        
        disJToI = sqrt((S(i).xd-S(PN_j).xd)^2 + (S(i).yd-S(PN_j).yd)^2 );
        if (disJToI < d0)
            Eij = ETX + Efs*(disJToI*disJToI);
        end
        if (disJToI >= d0)
            Eij = ETX + Emp*(disJToI*disJToI*disJToI*disJToI);
        end
        
        TEij = (Eij + ETX)*S(i).number_worker*bit;
        
        phi = (pi - (pi/2*S(PN_j).RE/S(PN_j).Initial_energy));
        Wij = exp(1/sin(phi));
        
        Cost_1 = Wij*TEij;
        
        for j= 2:1:length(PN)
            PN_j = PN(j);
            disJToI = sqrt((S(i).xd-S(PN_j).xd)^2 + (S(i).yd-S(PN_j).yd)^2 );
            if (disJToI < d0)
                Eij = ETX + Efs*(disJToI*disJToI);
            end
            if (disJToI >= d0)
                Eij = ETX + Emp*(disJToI*disJToI*disJToI*disJToI);
            end
 
            TEij = (Eij + ETX)*S(i).number_worker*bit;
 
            phi = (pi - (pi/2*S(PN_j).RE/S(PN_j).Initial_energy));
            Wij = exp(1/sin(phi));
 
            Costij = Wij*TEij;
            
            if (Costij < Cost_1) 
                parent = PN_j; 
                Cost_1 = Costij;
            end
        end
        
        S(i).parent = parent;
        S(i).Cost_parent = Cost_1;
    end
end
% --- End Routing
 %reduce energy
for i=1:1:n
    if (isequal(S(i).type,'W') && (S(i).RE >0))
        CH = S(i).CH;
        disToCH = sqrt((S(i).xd-S(CH).xd)^2 + (S(i).yd-S(CH).yd)^2 );
        if (disToCH < d0)
            S(i).RE = S(i).RE - bit*(ETX + Efs*(disToCH*disToCH));
        end
        if (disToCH >= d0)
            S(i).RE = S(i).RE - bit*(ETX + Emp*(disToCH*disToCH*disToCH*disToCH));
        end
    end
    
    if (isequal(S(i).type,'CH') && (S(i).RE >0))
        % NL nhan tu worker
        S(i).RE = S(i).RE - S(i).number_worker*bit*ETX;
        %Nl nhan tu CH khac
        parent = S(i).parent;
        % parent cua node i se bi mat di S(i).Cost_parent
        S(parent).RE = S(parent).RE - S(i).Cost_parent;
        % node i se mat them nl phat toi parent cua no
        disToParent = sqrt((S(i).xd-S(parent).xd)^2 + (S(i).yd-S(parent).yd)^2 );
        if (disToParent < d0)
            S(i).RE = S(i).RE - S(i).number_worker*bit*(ETX + Efs*(disToParent*disToParent));
        end
        if (disToParent >= d0) 
            S(i).RE = S(i).RE - S(i).number_worker*bit*(ETX + Emp*(disToParent*disToParent*disToParent*disToParent));
        end
    end
end
end
