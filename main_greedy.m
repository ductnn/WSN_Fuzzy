clear all;
close all;
clc;

% figure(1);

load('wsn.mat');
% plot(S(n+1).xd,S(n+1).yd,'x');

figure(1);
bit = 4000;

for r=1:1:10
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
           % fis1 = readfis('dis_Fuzzyfitness1');
           Energy_level = S(i).RE/S(i).Initial_energy;
           S(i).Fuzzy_fitness1 = evalfis([Energy_level S(i).distoBS], fis1);

           % fis2 = readfis('dis_Fuzzyfitness2');
           S(i).Fuzzy_fitness2 = evalfis([S(i).neigh_des S(i).neigh_cost], fis2);

           % random alpha from [0.9 1]
           % alpha = rand(1,1) / 10 + 0.9;
           alpha = 0.9423;
           
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
           % fis3 = readfis('Cluster.radius');
           S(i).rad = evalfis([S(i).Fuzzy_fitness1 S(i).Fuzzy_fitness2], fis3);

           S(i).number_worker = 0;
           cluster = cluster +1;

           plot(S(i).xd,S(i).yd,'k*');
           
           % plot(S(i), S(i-1), '-');
           % compute node j received from i
           for t= 1:1:n
             if ((isequal(S(t).type,'N') || isequal(S(t).type,'W'))&& (S(t).RE >0))
                disJToI = sqrt( (S(i).xd-S(t).xd)^2 + (S(i).yd-S(t).yd)^2 );
                if (disJToI <= S(i).rad)
                    k = length(S(t).candidate) + 1;
                    S(t).type = 'W';
                    S(t).candidate(k) = i;
                    
                    
                    % plot([S(i).xd,S(t).xd], [S(i).yd, S(t).yd], 'red');
                    
                end
                
             end  
             % plot([S(i).xd,S(t).xd], [S(i).yd, S(t).yd], 'red');
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
    
    %----Begin Routing----
    for i=1:1:n
        if S(i).type == 'CH';
            distCHtoBS = sqrt((sink.x - S(i).xd)^2 + (sink.y - S(i).yd)^2);
            plot([S(i).xd sink.x], [S(i).yd sink.y]);
            max_D = 0;
            for j=1:1:n
                if (isequal(S(j).type,'W') && (S(j).RE > 0))
                    u = [sink.x - S(i).xd, sink.y - S(i).yd];
                    v = [S(j).xd - S(i).xd, S(j).yd - S(i).yd];
                    CosTheta = dot(u,v)/(norm(u)*norm(v));
                    
                    distCHtoW = sqrt((S(i).xd - S(j).xd)^2 + (S(i).yd - S(j).yd)^2);
                    dist_shadow = distCHtoW*CosTheta;
                    
                    if (distCHtoW <= Rmax && CosTheta > cos(30*pi/180) && dist_shadow <= Rmax )
                       
                        if (max_D < dist_shadow)
                            max_D = dist_shadow;
                            plot([S(i).xd S(j).xd], [S(i).yd S(j).yd]);
                            plot([S(j).xd sink.x], [S(j).yd sink.y]);
                            %plot([S(i).xd S(j).xd sink.x], [S(i).yd S(j).yd sink.y]);
                        end
                    end
                end
            end
            % plot([S(i).xd, sink.x], [S(i).yd, sink.y]);
        end
    end
    %----End Routing------
    
    
end