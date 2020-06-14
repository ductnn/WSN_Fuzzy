
function [path,cost,heuristic,iterations] = greedybfs(source,target,weights,heuristics,names,startNode,goalNode)

% GREEDYBFS performs greedy best first search on graph with source, target, 
% weights and heuristics vectors.
% 
% Syntax:     
% 
% [path,cost,heuristic,iterations] = greedybfs(source,target,weights,heuristics,names,startNode,goalNode)
% [path,cost,heuristic,iterations] = greedybfs(source,target,weights,heuristics,startNode,goalNode)
% 
% Inputs:
% 
% source = Vector or cell array containing starting nodes of each of the edge.
% target = Vector or cell array containing ending nodes of each of the edge.
% weights = Vector containing weights of each of the edge.
% heuristics = Vector containing heuristic values for each node (usually straight line distances).
% names = Cell array containing string names of each of the node.
% startNode = Initial node in the graph.
% goalNode = Goal node in the graph.
% 
% Outputs:
% 
% path = Cell array containing search path.
% cost = Cost of path returned.
% heuristic = Heuristic value of last node.
% iterations = Table containing greedybfs iteration summary.
% 
% Example 01:
% 
% s = {'A','A','A','B','B','C'};
% t = {'B','C','D','E','F','G'};
% w = [1 5 3 4 5 9];
% h = [5 2 3 6 4 1 0];
% [path,cost,heuristic,iterations] = greedybfs(s,t,w,h,'A','G')
% 
% Example 02:
% 
% s = [1 1 1 2 2 3];
% t = [2 3 4 5 6 7];
% w = [1 5 3 4 5 9];
% h = [5 2 3 6 4 1 0];
% names = {'A','B','C','D','E','F','G'};
% [path,cost,heuristic,iterations] = greedybfs(s,t,w,h,n,'A','G')
% 
% Example 03:
% 
% s = [1 1 1 2 2 3];
% t = [2 3 4 5 6 7];
% w = [1 5 3 4 5 9];
% h = [5 2 3 6 4 1 0];
% [path,cost,heuristic,iterations] = greedybfs(s,t,w,h,1,7)
% 
% Coded by Ali Asghar Manjotho
% Lecturer, CSE-MUET
% Email: ali.manjotho.ali@gmail.com


    iterations = table;
    
    % Refactor source, target, names, startNode & goalNode vectors to numbers
    
    % If names argument missing 
    if (nargin<7)
        
        % Fifth argument (i.e. names) is starting Node
        ssNode = names;
        
        % Sixth argument (i.e. startNode) is goal Node
        ggNode = startNode;
        
        [s,t,n,sNode,gNode] = refactor(source,target,weights,ssNode,ggNode);
        
    else
        
        % Sixth argument (i.e. startNode) is starting Node
        ssNode = startNode;
        
        % Seventh argument (i.e. goalNode) is goal Node
        ggNode = goalNode;
        
        [s,t,n,sNode,gNode] = refactor(source,target,weights,names,ssNode,ggNode);        
        
    end

    
    % Get all unique nodes from source and target vectors
    uniqueNodes = getNodes(s,t);
    
    % Priority queue
    queue = [];
    
    % Initial path from starting node
    path = struct('Path',sNode,'Cost',0,'Heuristic',heuristics(sNode));
    
    % Add initial path to priority queue
    queue = [queue path];
    
    
    % Local variables to track iteration number
    iteration = 1; 
    
    % Update Iterations table
    iterations = [iterations; updateTable(s,t,n,queue,iteration)];
    
    
    
    % Repeat until queue is empty or goal is reached
    while(isGoalReached(queue,gNode)==0 && length(queue)>0)        
 
        % Put empty table row
        array = {'_____','_________________________________','_____','_____'};
        iterations = [iterations; cell2table(array)];
    
    
        % Get and remove minimum path from priority queue
        [minI,minP] = minPath(queue);
        queue(minI) = [];
        
        %Generate new paths
        newPaths = getNewPaths(s,t,weights,heuristics,minP);
        
        queue = [queue newPaths]; 
        
        % Update Iterations table
        iteration = iteration + 1; 
        iterations = [iterations; updateTable(s,t,n,queue,iteration)];
    
         
    end
    
    if(length(queue)>0)
    
        [minI,minP] = minPath(queue);
        path = n(minP.Path);
        cost = minP.Cost;
        heuristic = minP.Heuristic;
        
                
    else
        
        path = [];
        cost = [];
        heuristic = [];
        
    end
    
    iterations.Properties.VariableNames = {'Iteration' 'PriorityQueue' 'Cost' 'Heuristic'};

        
end



function [minIndex,path] = minPath(paths)

    minIndex = [];
    path = [];
    
    if(length(paths)>0)
        
        minIndex = 1;
        path = paths(minIndex);
    
        if(length(paths)>1)
           
            for i=2:length(paths)
                
                if(paths(i).Heuristic < path.Heuristic)
                
                    minIndex = i;
                    path = paths(minIndex);
                
                end
                
            end
            
        end
    
    end

end



function isGoal = isGoalReached(paths,goalNode)

    if(length(paths)==0)
        isGoal = 0;
        return;
    end
    
    [minI,minP] = minPath(paths);
    
    if(minP.Path(length(minP.Path)) == goalNode)
        isGoal = 1;
    else
        isGoal = 0;
    end

end


function weight = getWeight(s,t,weights,nodeA,nodeB)
    
    for i=1:length(s)
       
        if(s(i)==nodeA && t(i)==nodeB)
            weight = weights(i);
        end
        
    end

end


function paths = getNewPaths(s,t,w,h,path)

    paths = [];
    
    uniqueNodes = getNodes(s,t);
    
    if(~isempty(path))
        currentNode = path.Path(length(path.Path));

        childs = getChilds(s,t,currentNode);

        for i=1:length(childs)

            % If path is not a loop
            
            if(length(find(path.Path==childs(i)))==0)
                c = path.Cost + getWeight(s,t,w,currentNode,childs(i));
                heur = h(find(uniqueNodes==childs(i)));
                p = struct('Path',[path.Path childs(i)],'Cost',c,'Heuristic',heur);

                paths = [paths p];
            end

        end
    end

end


function childs = getChilds(source,target,node)
    
    childs = sort(target(find(source==node)));
    
end



function nodes = getNodes(s,t)

    nodes = unique(horzcat(s,t));

end


function [s,t,n,sn,gn] = refactor(source,target,weights,names,startNode,goalNode)

    % If names argument missing 
    if (nargin<6)
        
        % Fourth argument (i.e. names) is starting node
        sn = names;
        
        % Fifth argument (i.e. startNode) is goal node
        gn = startNode;

    else
        
        % Fifth argument (i.e. startNode) is starting node
        sn = startNode;  
        
        % Sixth argument (i.e. goalNode) is goal node
        gn = goalNode;
        
    end

    
    % Get all unique nodes
    uNodes = unique(horzcat(source,target));
        
        
    
    % If source and target are cell arrays
    if(iscell(source) && iscell(target))
    
        % If names argument missing
        if(nargin<6)
            n = uNodes;
        else
            n = names;
        end
        
        
        % Get unique nodes cell array
        uNodes = unique(horzcat(source,target));

        s = [];
        t = [];

        % Populate source and target with equivalent numeric values

        for i=1:length(source)
            [sFound,sIndex] = ismember(source(i),uNodes);
            [tFound,tIndex] = ismember(target(i),uNodes);
            s = [s sIndex];
            t = [t tIndex];
        end
            
        
        
        
    else
        
        s = source;
        t = target;
        
        % If names argument missing
        if(nargin<6)    
            
            uNodes = unique(horzcat(source,target));
            n = cell(1,length(uNodes));
            
            
            for i=1:length(uNodes)
                n{i} = num2str(uNodes(i));
            end
            
        else
            n = names;
        end
    end
    
    
    
    % If starting node is not a number
    if(~isnumeric(sn))

        sn = find(ismember(n,sn));
        
    end
    
    % If goal node is not a number
    if(~isnumeric(gn))

        gn = find(ismember(n,gn));
        
    end

end



function tableIteration = updateTable(s,t,n,queue,iteration)

    tempTable = table;

    uniqueNodes = getNodes(s,t);
    
    unsortedH = [];
    sortedH = [];



    for i=1:length(queue)
       unsortedH = [unsortedH queue(i).Heuristic]; 
    end

    mx = max(unsortedH);

    while(length(sortedH) ~= length(unsortedH))

        [mins,indices] = min(unsortedH);

        for j=1:length(indices)
            unsortedH(indices(j))=mx+1;
        end

        sortedH = [sortedH indices];

    end
    
    
    % Display current queue
    
    for p = 1:length(queue)
       
        path = queue(sortedH(p));
        
        pathStr = '<';

        for i=length(path.Path):-1:1
            if(i==1)
                pathStr = strcat(pathStr,sprintf('%s',char(n(find(uniqueNodes==path.Path(i))))));
            else
                pathStr = strcat(pathStr,sprintf('%s,',char(n(find(uniqueNodes==path.Path(i))))));
            end
        end
        pathStr = strcat(pathStr,sprintf('>'));



        % Display path cost    

        cost = path.Cost;
        
        % Display path heuristic    

        heuristic = path.Heuristic;
        
        array = {num2str(iteration) pathStr num2str(cost) num2str(heuristic)};
        tempTable = [tempTable; cell2table(array)];

    end
    
    tableIteration = tempTable;


end