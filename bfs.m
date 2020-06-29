function [searchNodes,iterations] = bfs(source,target,names,startNode)

%    BFS performs breadth first search on graph with source and target
%    vectors.
%     
%    Syntax:     
%          
%    [searchNodes,iterations] = bfs(source,target,names,startNode)
%    [searchNodes,iterations] = bfs(source,target,startNode)
%    
%    Inputs:
%
%    source = Vector or cell array containing starting node of each of the edge.
%    target = Vector or cell array containing ending node of each of the edge.
%    names = Cell array containing string names of each of the node.
%    startNode = Initial node in the graph.
%    
%    Outputs:
%
%    path = Cell array containing search path.
%    iterations = Table containing bfs iteration summary.
% 
%    Example 01:
%     
%    s = {'A','A','A','B','B','C'};
%    t = {'B','C','D','E','F','G'};
%    [searchNodes,iterations] = bfs(s,t,'A')
%     
%    Example 02:
%     
%    s = [1 1 1 2 2 3];
%    t = [2 3 4 5 6 7];
%    names = {'A','B','C','D','E','F','G'};
%    [searchNodes,iterations] = bfs(s,t,names,'A')
%     
%    Example 03:
%     
%    s = [1 1 1 2 2 3];
%    t = [2 3 4 5 6 7];
%    [searchNodes,iterations] = bfs(s,t,1)
%     
%    Coded by Ali Asghar Manjotho
%    Lecturer, CSE-MUET
%    Email: ali.manjotho.ali@gmail.com

    iterations = table;
    
    % Refactor source, target, names & startNode vectors to numbers
    
    % If names argument missing 
    if (nargin<4)
        
        % Third argument (i.e. names) is starting Node
        startingNode = names;
        
        [s,t,n,sNode] = refactor(source,target,startingNode);
        
    else
        
        % Fourth argument (i.e. startNode) is starting Node
        startingNode = startNode;
        
        [s,t,n,sNode] = refactor(source,target,names,startingNode);        
        
    end

    
    % Get all unique nodes from source and target vectors
    uniqueNodes = getNodes(s,t);
     
    % Initialize visited list and queue
    visited = [];
    queue = [];
     
    % Set starting node as current node and add it in to visited list
    currentNode = sNode;
    visited = [visited sNode];
          
    % Local variables to track iteration number
    iteration = 1;   

    
    % Update Iterations table
    iterations = [iterations; updateTable(s,t,n,currentNode,queue,visited,iteration,'Starting Node')];
    
    
    % Repeat until queue is empty
    while(~isempty(currentNode))
        
        % Get all childs of current node
        childs = getChilds(s,t,currentNode);
        
            
        for i=1:length(childs)
           
            

            % If new unvisited child found add it in queue and visited list
            if(length(find(visited==childs(i)))==0)

                queue = [queue childs(i)];
                visited = [visited childs(i)];
                
                % Increase iteration number
                iteration = iteration + 1;
                iterations = [iterations; updateTable(s,t,n,currentNode,queue,visited,iteration,strcat('Unvisited node found ',n(childs(i))))];

            end          

        end
         
              
         
        % If no new child found for current node then remove first item
        % from queue and make it as current node
        if (length(queue)>0)            
            currentNode = queue(1);
            queue(1) = [];
            
            
            comments = strcat('Dequeue ',n(currentNode),' from queue');
            
            
            iteration = iteration + 1;
            iterations = [iterations; updateTable(s,t,n,currentNode,queue,visited,iteration,comments)];
        else
            currentNode = [];
        end
        
               
         
    end
    
    
    % Update iteration table for last state of queue, visited list, current
    % node and search path when after queue is empty
    if(iteration > 0)
        
        iteration = iteration +1;
        
        iterations = [iterations; updateTable(s,t,n,currentNode,queue,visited,iteration,'BFS Converged')];

    end
    
    searchNodes = n(visited);
    iterations.Properties.VariableNames = {'Iteration' 'CurrentNode' 'Queue' 'Visited' 'Comments'};
    
        
end




function childs = getChilds(source,target,node)
    
    childs = sort(target(find(source==node)));
    
end



function nodes = getNodes(s,t)

    nodes = unique(horzcat(s,t));

end


function [s,t,n,sn] = refactor(source,target,names,startNode)

    % If names argument missing 
    if (nargin<4)
        
        % Third argument (i.e. names) is starting Node
        sn = names;

    else
        
        % Fourth argument (i.e. startNode) is starting Node
        sn = startNode;  
        
    end

    
    % Get all unique nodes
    uNodes = unique(horzcat(source,target));
        
        
    
    % If source and target are cell arrays
    if(iscell(source) && iscell(target))
    
        % If names argument missing
        if(nargin<4)
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
        if(nargin<4)    
            
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

end



function tableIteration = updateTable(s,t,n,currentNode,queue,visited,iteration,comments)

   
    uniqueNodes = getNodes(s,t);


    % Display current queue

    queueStr = '[';
    
    for i=length(queue):-1:1
        if(i==1)
            queueStr = strcat(queueStr,sprintf('%s',char(n(find(uniqueNodes==queue(i))))));
        else
            queueStr = strcat(queueStr,sprintf('%s ',char(n(find(uniqueNodes==queue(i))))));
        end
    end
    queueStr = strcat(queueStr,sprintf(']'));
    
    

    % Display current visited list     

    visitedStr = '[';
    
    for i=1:length(visited)
        if(i==length(visited))
            visitedStr = strcat(visitedStr,sprintf('%s',char(n(find(uniqueNodes==visited(i))))));
        else
            visitedStr = strcat(visitedStr,sprintf('%s ',char(n(find(uniqueNodes==visited(i))))));
        end
    end
    
    visitedStr = strcat(visitedStr,sprintf(']'));


    % Display current node
    if(~isempty(currentNode))
        node = n(currentNode);        
        currentNodeStr = sprintf('[%s]',node{1,1}(1,1));
    else
        currentNodeStr = sprintf('[]');
    end

    array = {iteration currentNodeStr queueStr visitedStr comments};
    tableIteration = cell2table(array);

end