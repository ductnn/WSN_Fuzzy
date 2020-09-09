clear all;
close all;
clc;

load('newtest.mat');
for r=1:1:4000
    r
    figure(1);
    hold on;
    
    dead = 0;
    for i=1:1:n
        %checking if there is a dead node
        if (S(i).RE<=0)
%             plot(S(i).xd,S(i).yd,'red D');
            S(i).RE = 0;
            dead = dead+1;
            S(i).state='DEAD';
            S(i).type = 'DEAD';
            S(i).number_worker = 0;
%             hold on;
        else
            S(i).type='N';
            S(i).state='Initial_State';
            S(i).number_worker = 0;
%             plot(S(i).xd,S(i).yd,'o');
%             hold on;
        end
    end
    S(n+1).id = n+1;
    S(n+1).xd = 50;
    S(n+1).yd = 50;
    S(n+1).Td = n+1;
%     plot(S(n+1).xd,S(n+1).yd,'x');
%     text(S(n+1).xd,S(n+1).yd,'  BS','Color','b','FontWeight','b');
    DEAD(r)=dead;
    for i= 1:1:n
        if S(i).RE > 0
            S(i).distoBS = norm([S(i).xd-S(n+1).xd S(i).yd-S(n+1).yd]);
            number_neighbor_i = 0;
            sigma_neigh_cost = 0;
            for j= 1:1:n
                disJToI = norm([S(i).xd-S(j).xd S(i).yd-S(j).yd]);
                if (disJToI <= Rmax && disJToI > 0)
                    number_neighbor_i = number_neighbor_i + 1;
                    sigma_neigh_cost = sigma_neigh_cost + disJToI^2;
                end
            end
            S(i).neigh_des  = number_neighbor_i/n;%HERE
            S(i).neigh_cost = (sqrt(sigma_neigh_cost/number_neighbor_i))/Rmax;
            %HERE
        end
    end
    
    %  Tinh Fuzzy va Td %
    for i= 1:1:n
        if (S(i).RE >0)
            % compute Td 
            Energy_level = S(i).RE/S(i).Initial_energy;
            S(i).Fuzzy_fitness1 = evalfis([Energy_level S(i).distoBS], fis1);
            S(i).Fuzzy_fitness2 = evalfis([S(i).neigh_des S(i).neigh_cost], fis2);
            alpha = rand(1,1)/10 + 0.9;
            S(i).Td = alpha * (1 - S(i).Fuzzy_fitness1) * Tc;       
            S(i).candidate = [];
%             plot(S(i).xd,S(i).yd,'o');
%             hold on;
        end
    end
    %Start bau CH
    %Sap xep S(i) theo chieu tang dan cua Td
    [x,idx] = sort([S.Td]);
    S = S(idx);
    for i=1:1:n
        if S(i).RE > 0
          if S(i).type == 'W'
              continue;
          else
              S(i).rad = evalfis([S(i).Fuzzy_fitness1 S(i).Fuzzy_fitness2],fis3);
              S(i).type = 'CH';
%               plot(S(i).xd,S(i).yd,'k*');
              for t=1:1:n
                distance = norm([S(i).xd-S(t).xd S(i).yd-S(t).yd]);
                if distance <= S(i).rad && distance > 0 && (S(i).RE>0)
                  k = length(S(t).candidate) + 1;
                  S(t).candidate(k) = S(i).id;
                  S(t).type = 'W';
                end
              end
          end
        end
    end
    % Remember to sort S
    [x,idx] = sort([S.id]);
    S = S(idx);
    for i=1:1:n
      if isequal(S(i).type,'W') && (S(i).RE >0) && ~isempty(S(i).candidate)
        CH_cost = 0;
        for t=1:1:length(S(i).candidate)
          distoCH = norm([S(i).xd-S(S(i).candidate(t)).xd S(i).yd-S(S(i).candidate(t)).yd]);
          CH_cost(t) = distoCH * S(S(i).candidate(t)).distoBS/S(S(i).candidate(t)).RE;
        end
        k = find(CH_cost==min(CH_cost));
        S(i).candidate = S(i).candidate(k);
      end
    end
    for i=1:1:n
        if isequal(S(i).type, 'CH') && S(i).RE >0
            x = [S.candidate];
            tf1 = x == S(i).id;
            S(i).number_worker = sum(tf1==1);
        end
    end
    %----End Cluster----
    
    %----Setup Routing-----
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
    neighbor = [];
    for i = 1:1:length(All_CH)
        neighbor = [neighbor;Func_Fuzzy_Greedy(All_CH,All_CH(i,1),30,fis)];
    end
    All_CH = [neighbor All_CH];
    for i = 1:1:length(All_CH)
        S(All_CH(i,2)).fuzzyrouting = All_CH(i,1);
    end
    for i=1:1:n
        if isequal(S(i).type,'W') && (S(i).RE >0)
                CH = S(i).candidate;
                distoCH = sqrt((S(i).xd-S(CH).xd)^2 + (S(i).yd-S(CH).yd)^2 );
                if distoCH < d0 
                   S(i).RE = S(i).RE - bit*(ETX + Efs*(distoCH^2));
                end
                if distoCH >= d0
                    S(i).RE = S(i).RE - bit*(ETX + Emp*(distoCH^4));
                end
        end
        if isequal(S(i).type,'CH') && (S(i).RE > 0)
            S(i).RE = S(i).RE - S(i).number_worker*bit*ETX;
            nextCH = S(i).fuzzyrouting;
            if nextCH ~= -1
                distEachCH = norm([S(i).xd-S(nextCH).xd S(i).yd-S(nextCH).yd]);
                if distEachCH < d0 
                    S(i).RE = S(i).RE - S(i).number_worker*bit*(ETX + Efs*(distEachCH^2));
                 end
                 if distEachCH >= d0
                     S(i).RE = S(i).RE - S(i).number_worker*bit*(ETX + Emp*(distEachCH^4));
                 end   
            else
                distoBS = norm([S(i).xd-S(n+1).xd S(i).yd-S(n+1).yd]);
                if distoBS < d0 
                    S(i).RE = S(i).RE - S(i).number_worker*bit*(ETX + Efs*(distoBS^2));
                 end
                 if distoBS >= d0
                     S(i).RE = S(i).RE - S(i).number_worker*bit*(ETX + Emp*(distoBS^4));
                 end
            end
        end
    end
    [S.RE]
end