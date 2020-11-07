clear All;
close all;
clc;
load('mazumdar2.mat');
% figure(1);
DEAD = 0;
S(n+1).yd = 50;
Rmax = 30;
for r= 1:1:2000
% Compute Neigbor desity & neighbor cost
    r
    dead = 0;
    for i=1:1:n
        %checking if there is a dead node
        if (S(i).RE <= 0)
%             plot(S(i).xd,S(i).yd,'red D');
            S(i).RE = 0;
            dead = dead+1;
            S(i).state = 'DEAD';
            S(i).type = 'DEAD';
            S(i).number_worker = 0;
%             hold on;
        else
            S(i).type = 'N';
            S(i).state = 'Initial_State';
            S(i).number_worker = 0;
            S(i).CH = 0;
%             plot(S(i).xd,S(i).yd,'o');
%             hold on;
        end
    end
    DEAD(r) = dead;
    for i=1:1:n
        if S(i).RE > 0
            S(i).distoBS = norm([S(i).xd-50 S(i).yd-50]);
            number_neighbor_i = 0;
            sigma_neigh_cost = 0;
            for j=1:1:n
                disJToI = norm([S(i).xd-S(j).xd S(i).yd-S(j).yd]);
                if (disJToI <= Rmax && disJToI > 0)
                    number_neighbor_i = number_neighbor_i + 1;
                    sigma_neigh_cost = sigma_neigh_cost + disJToI^2;
                end
            end
            S(i).neigh_des  = number_neighbor_i/n;
            S(i).neigh_cost = (sqrt(sigma_neigh_cost/number_neighbor_i))/Rmax;
        end
    end
    
    %  Tinh Fuzzy va Td %
    for i= 1:1:n
        if (S(i).RE >0 )
            % compute Td 
            Energy_level = S(i).RE/S(i).Initial_energy;
            S(i).Fuzzy_fitness1 = evalfis([Energy_level S(i).distoBS], fis1);
            % S(i).Fuzzy_fitness2 = evalfis([S(i).neigh_des S(i).neigh_cost], fis2);
            alpha = rand(1,1)/10 + 0.9;
            S(i).Td = alpha * (1 - S(i).Fuzzy_fitness1) * Tc;       
            S(i).candidate = [];
%             plot(S(i).xd,S(i).yd,'o');
%             hold on;
        end
    end
    
    %Start bau CH
    %Sap xep S(i) theo chieu tang dan cua Td

    [x, idx] = sort([S.Td]);
    S = S(idx);
    sink.level = [15 20 25 30];
    for i=1:1:n
        if S(i).RE > 0
          if isequal(S(i).type, 'W')
              continue;
          else
            if isequal(S(i).type, 'N')
              for t=1:1:length(sink.level)
                    if S(i).distoBS <= 10 
                        S(i).rad = 0.9 * 13;
                    elseif S(i).distoBS <= 15 && S(i).distoBS > 10
                        S(i).rad = 0.9 * 8;
                    elseif S(i).distoBS <= 20 && S(i).distoBS > 15
                        S(i).rad = 0.9 * 5; 
                    elseif S(i).distoBS <= 25 && S(i).distoBS > 20
                        S(i).rad = 0.9 * 3;
                    elseif S(i).distoBS <= 30 && S(i).distoBS > 25
                        S(i).rad = 0.9 * 2;
                    else
                        S(i).rad = 0.9 * 1;
                    end
              end
            end
            %   S(i).rad = evalfis([S(i).Fuzzy_fitness1 S(i).Fuzzy_fitness2],fis3);
              S(i).type = 'CH';
%               plot(S(i).xd,S(i).yd,'k*');
              for t=1:1:n
                distance = norm([S(i).xd-S(t).xd S(i).yd-S(t).yd]);
                if distance <= S(i).rad && distance > 0 && S(i).RE>0 && ~isequal(S(t).type, 'CH')
                    S(t).candidate = [S(t).candidate S(i).id];
                    S(t).type = 'W';
                end
              end
          end
        end
    end

    %Remember to sort S
    % CheckKKKKKKKKKKKKKKKKs
    hamlon = 0;
    [x,idx] = sort([S.id]);
    S = S(idx);
    for i=1:1:n
        if isequal(S(i).type, "W")
            for z=1:1:length(S(i).candidate)
                list_RE = [S(S(i).candidate).RE];
                [result, index] = max(list_RE);
                S(i).candidate = S(S(i).candidate(index)).id;
                S(i).CH = S(i).candidate;
            end
        end
    end
    % %==========================FIBONACCI=============================================
    % % Set sink level
    % % R_priv: Bán kính riêng của từng node


    % sink.level = 3.85 * [1 2 3 5 8 13];
    % % sink.level = [4 7.6 11.55 19.2 30.7 50.05];
    % for i=1:1:n
    %     if  isequal(S(i).type, 'CH')
    %         for t=1:1:length(sink.level)
    %             if S(i).distoBS <= sink.level(t)
    %                 S(i).R_priv = sink.level(t);
    %             end
    %         end

    %         % Append CH to candidate
    %         for j=1:1:n
    %             if S(j).type == 'W'
    %                 distTOCH = norm([S(i).xd-S(j).xd S(i).yd-S(j).yd]);
    %                 if distTOCH < S(i).R_priv
    %                     S(j).candidate = [S(j).candidate S(i).id];
    %                 end
    %             end
    %         end

    %         % Chon CH theo nang luong
    %         for k=1:1:n
    %             if S(k).type == 'W'
    %                 for z=1:1:length(S(k).candidate) 
    %                     x = [S(S(k).candidate).RE];
    %                     [a b] = max(x);
    %                     S(k).candidate = S(S(k).candidate(b)).id;
    %                     S(k).CH = S(k).candidate;
    %                 end
    %             end
    %         end

    %     end
    % end

    %==========================END_FIBONACCI=========================================

    % Add worker's number of each CH
    for i=1:1:n
        number_worker = 0;
        if isequal(S(i).type,'CH') && S(i).RE > 0
            for j=1:1:n
                if isequal(S(j).type, 'W') && S(i).RE > 0 
                    if S(j).candidate == S(i).id
                        number_worker = number_worker + 1;
                    end
                end
            end
        end
        S(i).number_worker = number_worker;
    end

    % --- End Bau CH

    %Algorithm 2: DFCR routing%

    % Calculate level L
    for i=1:1:n
        S(i).PN = [];
        if isequal(S(i).type,'CH') && S(i).RE > 0
            S(i).L = ceil(S(i).distoBS/Rmax);
        end
    end

    for i=1:1:n
        Cost_election = [];
        if isequal(S(i).type,'CH') && S(i).RE > 0
            for k=2:1:S(i).L
                for j=1:1:n
                    if isequal(S(j).type,'CH') && S(j).RE > 0
                        distance = norm([S(i).xd-S(j).yd S(i).yd-S(j).yd]);
                        if S(i).L > S(j).L && distance <= k * Rmax && S(i).id ~= S(j).id && distance > 0
                            if ~ismember(S(j).id, S(i).PN)

                                S(i).PN(length(S(i).PN)+1) = S(j).id;
                                if distance < d0
                                    Eij = ETX + Efs * (distance^2);
                                else
                                    Eij = ETX + Emp * (distance^4);
                                end

                                TEij = (Eij + ETX)*S(i).number_worker * bit;
                                phi = (pi - (pi * S(j).RE / (2 * S(j).Initial_energy)));
                                Wij = exp(1 / (sin(phi)));
                                Costij = TEij * Wij;
                                Cost_election = [Cost_election;Costij j TEij];
                                
                            end
                        end
                    end
                end
            end
            if ~isempty(Cost_election)
                [a,ida] = min(Cost_election(:,1));
                S(i).PN = S(Cost_election(ida,2)).id;
                S(i).cost = a;
                S(i).TEij = Cost_election(ida,3);
            end
        end
    end
    % Routing
    %End of routing
    for i=1:1:n
        if isequal(S(i).type,'W') && S(i).RE > 0

            CH = S(i).CH;
            distoCH = norm([S(i).xd-S(CH).xd S(i).yd-S(CH).yd]);
            if distoCH < d0
                S(i).RE = S(i).RE - bit * (ETX + Efs * (distoCH^2));
            else
                S(i).RE = S(i).RE - bit * (ETX + Emp * (distoCH^4));
            end

        elseif isequal(S(i).type,'CH') && S(i).RE > 0
            %Nang luong nhan
            S(i).RE = S(i).RE - (S(i).number_worker + 1) * bit * ETX;
            parent = S(i).PN;
            %Nang luong truyen
            if ~isempty(parent)
                distance = norm([S(i).xd-S(parent).xd  S(i).yd-S(parent).yd]);
%                 S(parent).RE = S(parent).RE - S(i).TEij * 10^-1;
                if distance < d0
                    S(i).RE = S(i).RE - (S(i).number_worker + 1) * bit *(ETX + Efs * (distance^2));
                else
                    S(i).RE = S(i).RE - (S(i).number_worker + 1) * bit * (ETX + Emp * (distance^4));                    
                end

            else

                if S(i).distoBS < d0
                    S(i).RE = S(i).RE - (S(i).number_worker + 1) * bit * (ETX + Efs * (S(i).distoBS^2));
                else
                    S(i).RE = S(i).RE - (S(i).number_worker + 1) * bit * (ETX + Emp * (S(i).distoBS^4));                    
                end
            end
        end
    end
    [S.RE]
end


