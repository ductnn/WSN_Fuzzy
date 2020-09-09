function neighbor = Func_Fuzzy_Greedy(All_CH,start_id,Rmax,fis4)
    neighbor = [];
    finish_point = [50 50];
    start_point = [All_CH(All_CH(:,1) == start_id,2) All_CH(All_CH(:,1) == start_id,3)];
    vector_sf = [finish_point(1)-start_point(1) finish_point(2)-start_point(2)];
    candidate_point = 0;
    if norm(vector_sf) >= Rmax
        for i=1:1:length(All_CH)-1
            vector_s_CH = [All_CH(i,2)-start_point(1) All_CH(i,3)-start_point(2)];
            %Tinh cos Theta cua 2 vector
            Costheta = (dot(vector_sf,vector_s_CH))/(norm(vector_s_CH)*norm(vector_sf));
            theta = acosd(Costheta);
            if (0 <= Costheta) && (Costheta <= 90) && 0 < norm(vector_s_CH) && norm(vector_s_CH)<= Rmax
                % Luu lai gia tri cua vector
                neighbor = [neighbor;All_CH(i,1) theta Costheta*norm(vector_s_CH)];
            end
        end
        if ~isempty(neighbor)
            Fuzzy_Routing = evalfis(neighbor(:,2:3),fis4);
            Fuzzy_Routing = [neighbor(:,1) Fuzzy_Routing];
            [a,b] = max(Fuzzy_Routing(:,2));
            candidate_point = [Fuzzy_Routing(b,1)];
            neighbor = candidate_point;
        else
            neighbor = -1;
        end
    else
        neighbor = 101;
    end
end