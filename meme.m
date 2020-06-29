clear all;
close all;
clc;

% figure(1);

load('test.mat');
Rmax = 30;
figure(1);

for r=1:1:9999
    r
    figure(1);
    hold off;
    
    dead = 0;
    for i=1:1:n
        %checking if there is a dead node
        if (S(i).RE<=0)
            plot(S(i).xd,S(i).yd,'red +');
            dead=dead+1;
            S(i).state='DEAD';
            S(i).type='DEAD';
            S(i).RE = 0;
            hold on;
        else
            S(i).type='N';
            S(i).state='Initial_State';
            plot(S(i).xd,S(i).yd,'o');
            hold on;
        end
    end
    plot(S(n+1).xd,S(n+1).yd,'x');
    text(S(n+1).xd,S(n+1).yd,'  BS','Color','b','FontWeight','b');
    
    STATISTICS(r+1).DEAD=dead;
    DEAD(r+1)=dead;
    
%     if(dead==1)
%         if(flag_first_dead==0)
%             first_dead=r;
%             flag_first_dead=1;
%             beep;
%             msgbox('Danger, A node has dead');
%         end
%     end
    
    %  CH %
    for i= 1:1:n
       if (S(i).RE >0 )
           % Compute Neigbor desity & neighbor cost
           number_neighbor_i = 0;
           sigma_neigh_cost = 0;
           for j= 1:1:n
               disJtoI = sqrt((S(i).xd-S(j).xd)^2 + (S(i).yd-S(j).yd)^2);
               if (disJtoI <= Rmax && disJtoI > 0)
                   number_neighbor_i = number_neighbor_i + 1;
                   sigma_neigh_cost = sigma_neigh_cost + disJtoI^2;
               end
           end
           S(i).neigh_des = number_neighbor_i/n; 
           S(i).neigh_cost = (sqrt(sigma_neigh_cost/number_neighbor_i))/Rmax;
           
           
           S(i).distoBS = sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2);

           % compute Td 
           Energy_level = S(i).RE/S(i).Initial_energy;
%            fis1 = readfis('dis_Fuzzyfitness1');           
           S(i).Fuzzy_fitness1 = evalfis([Energy_level S(i).distoBS], fis1);
           S(i).Fuzzy_fitness2 = evalfis([S(i).neigh_des S(i).neigh_cost], fis2);

           % random alpha from [0.9 1]
           alpha = rand(1,1) / 10 + 0.9;
%            alpha = 0.9423;
           
           S(i).Td = alpha * (1 - S(i).Fuzzy_fitness1) * Tc;
           %%S(i).rad = evalfis([S(i).RE S(i).distoBS S(i).Fuzzy_fitness1], fis1);

           %Initial candidates
           S(i).candidate = [];
       end
    end
    
    %disp([S.Td]);

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
               S(i).rad = evalfis([S(i).Fuzzy_fitness1 S(i).Fuzzy_fitness2], fis3);

               S(i).number_worker = 0;
               cluster = cluster +1;

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
    %----End Cluster----
    
    %----Begin Routing-----
    % countCH = 0;
    CH_number = S(strcmp({S.type},'CH'));
    %Hoanh do cua CH
    x_CH = [CH_number.xd,S(n+1).xd]';
    %Tung do cua CH
    y_CH = [CH_number.yd,S(n+1).yd]';
    %id cua CH (sap xep theo so thu tu cua ma tran input cua code)
    id_CH = [CH_number.id, 101]';
    %Trong nay chua toan bo diem CH va diem sink ben duoi cung
    All_CH = [id_CH x_CH y_CH];
    %Khai bao so luong cua cac node ket noi trong network
    num_nodes = length(All_CH);
    num_segs=0;
    
    segments = zeros(length(All_CH)*(length(All_CH)-1)/2,3);
    for a = 1:length(All_CH)-1 % create edges between some of the nodes
        text(All_CH(a,2),All_CH(a,3),[' ' num2str(id_CH(a))],'Color','b','FontWeight','b')
        for b = a+1:length(All_CH)
            d = sqrt(sum((All_CH(a,2:3) - All_CH(b,2:3)).^2));
            if (d < Rmax)
                %plot([All_CH(a,2) All_CH(b,2)],[All_CH(a,3) All_CH(b,3)],'k.-')
                % add this link to the segments list
                num_segs = num_segs + 1;
                segments(num_segs,:) = [num_segs All_CH(a,1) All_CH(b,1)];
            end

        end
    end
    %Like a routing table
    segments(num_segs+1:num_nodes * (num_nodes-1)/2,:) = [];
    %------End Routing--------
    
    %Reduce energy
    %Initial Energy bit
    Eb = 13e-9;
%      Eb=1e-6;
%     Energy_Transmission = 0;
    for i=1:1:n
        if (isequal(S(i).type,'W') && (S(i).RE >0))
            CH = S(i).CH;
            disToCH = sqrt((S(i).xd-S(CH).xd)^2 + (S(i).yd-S(CH).yd)^2 );
            if (disToCH < d0)
                S(i).RE = S(i).RE - bit*(ETX + Efs*(disToCH^2));
            end
            if (disToCH >= d0)
                S(i).RE = S(i).RE - bit*(ETX + Emp*(disToCH^4));
            end
        end
    end
    for i = 1:1:length(CH_number)
        [distance, path] = dijkstra(All_CH, segments, CH_number(i).id, 101);
%         path = Greedy(All_CH,CH_number(i).id,30,Rmax);
        Energy_Transmission = 0;
        
        if isnan(path)
            distanceNodetoCH = norm([CH_number(i).xd-50 CH_number(i).yd-50]);            
%             S(i).RE = S(i).RE -( (ETX+EDA)*bit + Efs*bit*(distanceNodetoCH^2 ));
            CH_number(i).RE = CH_number(i).RE - CH_number(i).number_worker*bit*(ETX+Efs*distanceNodetoCH^2);
            continue;
            
%             if (distanceNodetoCH >= d0)
%                 S(i).RE = S(i).RE -( (ETX+EDA)*bit + Efs*bit*(distanceNodetoCH^4 ));
%             end
            
        end
        for k = 1:1:length(path)-1
            Energy_Transmission = CH_number([CH_number.id] == path(k)).number_worker*Eb*bit + Energy_Transmission;
            CH_number([CH_number.id] == path(k)).RE = CH_number([CH_number.id] == path(k)).RE - Energy_Transmission;
            
            if(CH_number([CH_number.id] == path(k)).RE <= 0)
               CH_number([CH_number.id] == path(k)).RE = 0;
               S([S.id] == CH_number([CH_number.id] == path(k)).id).state = 'DEAD';
               S([S.id] == CH_number([CH_number.id] == path(k)).id).type = 'DEAD';
            end
            S([S.id] == CH_number([CH_number.id] == path(k)).id).RE = CH_number([CH_number.id] == path(k)).RE;
            
            
        end
    end
    % S.RE = CH_number.RE
%     [CH_number.RE]
end