clear All;
close all;
clc;
load('fibo_test_1.mat');
% figure(1);
DEAD = 0;
S(n+1).yd = 50;
Rmax = 30;
for r= 1:1:2000
    r
    dead = 0;
    for i=1:1:n
        %Kiểm tra nếu có nút chết
        if (S(i).RE <= 0)
            % plot(S(i).xd,S(i).yd,'red D');
            S(i).RE = 0;
            dead = dead+1;
            S(i).state = 'DEAD';
            S(i).type = 'DEAD';
            S(i).number_worker = 0;
            % hold on;
        else
            % Các nút quay lại trạng thái Normal node
            S(i).type = 'N';
            S(i).state = 'Initial_State';
            S(i).number_worker = 0;
            S(i).CH = 0;
            % plot(S(i).xd,S(i).yd,'o');
            % hold on;
        end
    end
    DEAD(r) = dead;
    %Tính toán neighbor density và neigh bor cost
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
            % Neighbor density và neighbor cost
            S(i).neigh_des  = number_neighbor_i/n;
            S(i).neigh_cost = (sqrt(sigma_neigh_cost/number_neighbor_i))/Rmax;
        end
    end
    
    %  Tính bộ Fuzzy Fitness 1 va Td %
    for i= 1:1:n
        if (S(i).RE >0 )
            % compute Td 
            Energy_level = S(i).RE/S(i).Initial_energy;
            S(i).Fuzzy_fitness1 = evalfis([Energy_level S(i).distoBS], fis1);
            % S(i).Fuzzy_fitness2 = evalfis([S(i).neigh_des S(i).neigh_cost], fis2);
            alpha = rand(1,1)/10 + 0.9;
            S(i).Td = alpha * (1 - S(i).Fuzzy_fitness1) * Tc;       
            S(i).candidate = [];
            % plot(S(i).xd,S(i).yd,'o');
            % hold on;
        end
    end
    
    %==============================Start bau CH==================================

    [x, idx] = sort([S.Td]);
    S = S(idx);
    % Bán kính của các nút CH sử dụng để chọn các Worker
    sink.level = 1.4 * [1 2 3 5 8 13 21];
    sink.radius = [5 10 15 20 25 30];
    %Chọn CH và Worker
    for i=1:1:n
        if S(i).RE > 0
          if isequal(S(i).type, 'W')
              continue;
          else
            if isequal(S(i).type, 'N')
                % Lựa chọn bán kính phù hợp cho mỗi vùng phủ bán kính từ điểm sink
              for t=1:1:length(sink.level)
                    if S(i).distoBS <= sink.radius(1) 
                        S(i).rad = sink.level(1);

                    elseif S(i).distoBS <= sink.radius(2) && S(i).distoBS > sink.radius(1)
                        S(i).rad = sink.level(2);

                    elseif S(i).distoBS <= sink.radius(3) && S(i).distoBS > sink.radius(2)
                        S(i).rad = sink.level(3); 

                    elseif S(i).distoBS <= sink.radius(4) && S(i).distoBS > sink.radius(3)
                        S(i).rad = sink.level(4);

                    elseif S(i).distoBS <= sink.radius(5) && S(i).distoBS > sink.radius(4)
                        S(i).rad = sink.level(5);

                    elseif S(i).distoBS <= sink.radius(6) && S(i).distoBS > sink.radius(5)
                        S(i).rad = sink.level(6);
                    elseif S(i).distoBS > sink.radius(6)
                        S(i).rad = sink.level(6);
                    end
              end
            end
            %   S(i).rad = evalfis([S(i).Fuzzy_fitness1 S(i).Fuzzy_fitness2],fis3);
              S(i).type = 'CH';
            %   plot(S(i).xd,S(i).yd,'k*');
              % Tìm Worker nằm trong bán kính của CH
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
    [x,idx] = sort([S.id]);
    S = S(idx);
    %Tìm nút CH có năng lượng nhiều nhất trong candidate của Worker
    for i=1:1:n
        if isequal(S(i).type, "W")
            for z=1:1:length(S(i).candidate)
                list_RE = [S(S(i).candidate).RE];
                % Tìm nút CH có năng lượng nhiều nhất của candidate
                [result, index] = max(list_RE);
                S(i).CH = S(S(i).candidate(index)).id;
            end
        end
    end
    %Tính số lượng Worker của mỗi CH
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

    % =======================End Bau CH===========================

    %Algorithm 2: DFCR routing%

    % Tính toán level L
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
                                %Tính Eij nhằm tính TEij
                                if distance < d0
                                    Eij = ETX + Efs * (distance^2);
                                else
                                    Eij = ETX + Emp * (distance^4);
                                end
                                %Tính TEij
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
    %===============Xong phần định tuyến===============

    %===============Tính toán năng lượng===============
    for i=1:1:n
        % Tính năng lượng truyền của các worker tới CH
        if isequal(S(i).type,'W') && S(i).RE > 0

            CH = S(i).CH;
            distoCH = norm([S(i).xd-S(CH).xd S(i).yd-S(CH).yd]);
            if distoCH < d0
                S(i).RE = S(i).RE - bit * (ETX + Efs * (distoCH^2));
            else
                S(i).RE = S(i).RE - bit * (ETX + Emp * (distoCH^4));
            end
        % Tính năng lượng truyền và năng lượng nhận của các nút CH
        elseif isequal(S(i).type,'CH') && S(i).RE > 0
            %Năng lượng nhận
            S(i).RE = S(i).RE - (S(i).number_worker + 1) * bit * ETX;
            parent = S(i).PN;
            %Năng lượng truyền trong trường hợp CH không có Parent Node
            if ~isempty(parent)
                distance = norm([S(i).xd-S(parent).xd  S(i).yd-S(parent).yd]);
%                 S(parent).RE = S(parent).RE - S(i).TEij * 10^-1;
            %Năng lượng truyền nếu khoảng cách tới Parent node < d0
                if distance < d0
                    S(i).RE = S(i).RE - (S(i).number_worker + 1) * bit *(ETX + Efs * (distance^2));
            %Năng lượng truyền nếu khoảng cách tới Parent node > d0
                else
                    S(i).RE = S(i).RE - (S(i).number_worker + 1) * bit * (ETX + Emp * (distance^4));                    
                end
            %Năng lượng truyền trong trường hợp CH có Parent Node
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


