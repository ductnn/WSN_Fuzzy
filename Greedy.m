function path = Greedy(All_CH,start_id,angle,Rmax)
    %Greedy
    %vector_finish = [xd yd]
    n = 1;
    %Toa do diem finish: [xd yd]
    distance_max = 0;
    finish_point = [50 50];
    path(n) = start_id;
    while true
        k = find(All_CH(:,1) == start_id);
        %Toa do diem start:[xd yd]
        start_point = [All_CH(find(All_CH(:,1) == start_id),2) All_CH(find(All_CH(:,1) == start_id),3)];
        %Vector cua diem start voi diem sink[x_sink-x_start]
        vector_sf = [finish_point(1)-start_point(1) finish_point(2)-start_point(2)];
        %Kiem Tra goc va do dai lon nhat cua cac diem xung quanh no
        n
        n = n + 1;
        if norm(vector_sf) <= Rmax
            path(n) = 101;
            break;
        end
        for i=1:1:length(All_CH)
            %Vector cua diem start voi tung diem trong All_CH
            vector_s_CH = [All_CH(i,2)-start_point(1) All_CH(i,3)-start_point(2)];
            %Tinh cos Theta cua 2 vector
            Costheta = (dot(vector_sf,vector_s_CH))/(norm(vector_s_CH)*norm(vector_sf));
            if Costheta >= cos(angle * pi/180) && norm(vector_s_CH) <= Rmax
                if distance_max < Costheta * norm(vector_s_CH)
                    distance_max = Costheta * norm(vector_s_CH);
                    %Luu gia tri id cua diem cao nhat
                    k = i;
                end

            end
        end
        distance_max=0;
        start_id = All_CH(k,1);
        if isempty(find(path==start_id, 1))
            path(n) = start_id;
        else
            break;
        end

    end
    if path(length(path)) ~= 101
        path = nan;
    end
end