classdef dllist < handle
    
    properties
        Head = dlnode.empty;
        Num_nodes = 0;
    end
    
    methods
        
        % construct an empty linked-list
        function list = dllist(varargin)
        end
       
        % adds a node to the list and updates the number of nodes
        function add_node(list, new_node)
            if (list.Num_nodes == 0)
                list.Head = new_node;
            else
                new_node.insertAfter(list.get_node(list.Num_nodes));
            end
            list.Num_nodes = list.Num_nodes + 1;
        end
        
        % removes a specified node (referenced by sequential number) from
        % the list and updates the number of nodes
        function remove_node(list, node_number)
            switch list.Num_nodes
                case 0
                    return;
                case 1
                    list.Head = dlnode.empty;
                    list.Num_nodes = 0;
                otherwise
                    switch node_number
                        case 1
                            list.Head = list.get_node(2);
                            list.Head.Prev.removeNode();
                        otherwise
                            node_to_remove = list.get_node(node_number);
                            node_to_remove.removeNode()
                    end
                    list.Num_nodes = list.Num_nodes - 1;
            end
        end
        
        % get the specified node (referenced by sequential number) in the
        % linked list
        function node = get_node(list, node_num)
            if (node_num <= list.Num_nodes)
                counter = 1;
                curr_node = list.Head;
                while counter < node_num
                    curr_node = curr_node.Next;
                    counter = counter + 1;
                end
                node = curr_node;
            else
                node = dlnode.empty;
            end
        end
        
        % clear the entire list of nodes following the dlnode.clearList
        % function
        function clearList(list)
            if list.Num_nodes == 0
                list.Head = dlnode.empty;
                list.Num_nodes = 0;
            else
                lastNode = get_node(list, list.Num_nodes);
                lastNode.clearList();
                list.Head = dlnode.empty;
                list.Num_nodes = 0;
            end
            
        end
        
    end
    
end

