classdef dlnode < handle

    properties (SetAccess = private)
        Next = dlnode.empty;
        Prev = dlnode.empty;
   end
    
    methods
        
        % Creating an individual node (not connected)
        function node = dlnode(varargin)
        end
       
        % Insert node into a doubly linked list after specified node, or
        % link the two specified nodes if there is not already a list.
        % Assigns the correct values for Next and Prev properties.
        function insertAfter(newNode, nodeBefore)
            removeNode(newNode);
            newNode.Next = nodeBefore.Next;
            newNode.Prev = nodeBefore;
            if ~isempty(nodeBefore.Next)
                nodeBefore.Next.Prev = newNode;
            end
            nodeBefore.Next = newNode;
        end
        
        % Insert node into doubly linked list before specified node, or
        % link the two specified nodes if there is not already a list. This
        % method assigns correct values for Next and Prev properties.
        function insertBefore(newNode, nodeAfter)
            removeNode(newNode);
            newNode.Next = nodeAfter;
            newNode.Prev = nodeAfter.Prev;
            if ~isempty(nodeAfter.Prev)
                nodeAfter.Prev.Next = newNode;
            end
            nodeAfter.Prev = newNode;
        end
        
        % Remove node and fix the list so that remaining nodes are properly
        % connected. node argument must be scalar. Once there are no
        % references to node, MATLAB deletes it.
        function removeNode(node)
            if ~isscalar(node)
                error('Nodes must be scalar')
            end
            prevNode = node.Prev;
            nextNode = node.Next;
            if ~isempty(prevNode)
                prevNode.Next = nextNode;
            end
            if ~isempty(nextNode)
                nextNode.Prev = prevNode;
            end
            node.Next = dlnode.empty;
            node.Prev = dlnode.empty;
        end
        
        % Avoid recursive calls to destructor as a result of clearing the
        % list variable. Loop through list to disconnect each node. When
        % there are no references to a node, MATLAB calls the class
        % destructor (see the delete method) before deleting it.
        function clearList(node)
            prev = node.Prev;
            next = node.Next;
            removeNode(node)
            while ~isempty(next)
                node = next;
                next = node.Next;
                removeNode(node);
            end
            while ~isempty(prev)
                node = prev;
                prev = node.Prev;
                removeNode(node)
            end
        end
        
    end
    
    methods (Access = private)
        
        % Class destructor method. MATLAB calls the delete method you
        % delete a node that is still connected to the list.
        function delete(node)
            clearList(node)
        end
        
    end
    
end

