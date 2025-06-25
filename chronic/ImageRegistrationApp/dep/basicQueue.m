classdef basicQueue < handle

    %THIS IS A VERY BASIC CIRCULAR QUEUE WHICH IS NOT OPTIMISED
    %EACH TIME THE QUEUE IS SHORTENED OR LENGHTENED THE ARRAY GETS
    %REWRITTEN FULLY
    %COULD LOOK AT A SOLUTION WHERE A MAX LENGTH IS ADDED, I.e. 64, or 256
    %etc... AND THEN KEEP AN INDEX OF THE END AND LENGTH OF THE "USABLE"
    %LIST, THEREFORE IF WE DYNAMICALLY SHORTEN AND LENGTHEN IT DOESN'T LOSE
    %TRACK OF WHERE THE "USABLE" LIST SHOULD LOOP BACK TO, ESPECIALLY WHEN
    %NOT YOU JUST ASSIGN AND IT LOOPS BACK AROUND THE QUEUE
    properties (Access = public)
        vals = [];
        back = 0; %back of queue, oldest vals
        size = 0; %total arr length
        front = 0; %front of queue, newest vals
        len = 0; %number of assigned vals
    end 
    methods (Access = private)
        
        function rewrite(obj, size)
            arguments
                obj
                size = -1
            end
            if size == -1
                size = obj.size;
            end

            tempVals = zeros(1, size);
            tLen = obj.len;
            
            for iter = 1:tLen
                tempVals(iter) = obj.pop();
            end

            obj.vals = tempVals;
            obj.back = 1;
            obj.front = tLen;
            obj.len = tLen;
            obj.size = size;
        end

        function add(obj, val)
            obj.front = mod(obj.front, obj.size) + 1;
            obj.vals(obj.front) = val;
            obj.len = obj.len + 1;
        end

        function val = isempty(obj)
            if obj.len == 0
                val = 1;
            else
                val = 0;
            end
        end
    end

    methods (Access = public)
        
        function bulkRemoveAndReassign(obj, selection)
            %exclusively used when imRegApp deletes rows, as such the ID's
            %get changed so I need to change it in this object
            selection = sort(selection);
            tselection = intersect(selection, obj.vals);
            
            %remove the rows with indices in the queue
            obj.rewrite()

            for i=1:length(tselection)
                idx = find(obj.vals == tselection(i));
                obj.vals(idx) = [];
                obj.len = obj.len - 1;
                obj.front = obj.front - 1;
                obj.size = obj.size - 1;
            end
            
            %shift indice value by the required amount
            for i=1:length(obj.vals)
                correction = sum(selection < obj.vals(i));
                obj.vals(i) = obj.vals(i) - correction;
            end

            if obj.size < obj.len
                disp(obj.size)
                disp(obj.len)
                disp(obj.getState)
                error("this is a bug with bulkRemoveAndReassign")
            end
        end

        function loadState(obj, state)
            obj.vals = state.vals;
            obj.back = state.back ;
            obj.size = state.size;
            obj.front = state.front;
            obj.len = state.len;
            obj.rewrite();
        end

        function state = getState(obj)
            state.vals = obj.vals;
            state.back = obj.back;
            state.size = obj.size;
            state.front = obj.front;
            state.len = obj.len;
        end

        function retVal = enqueue(obj, val)
            retVal = 0;
            
            if obj.size == 0
                error("cannot enqueue as queue size is 0, please use lengthen(<array of integers to add>) or setlen(<length of queue>)")
            end

            if obj.len == obj.size
                retVal = obj.pop();
            end

            obj.add(val);
        end

        function lengthen(obj, val)
            %increase length of queue by adding giving an array of vals to
            %add to the head
            arguments
                obj 
                val double = 0 %array of vals 
            end
            obj.rewrite(obj.size + length(val));

            for i=1:length(val)
                obj.add(val(i));
            end
        end

        function retVal = setLen(obj, val)
            %increase length of queue by specifying new max length
            retVal = 0;
            
            if val < 0
                error("length must be 0 or greater than 0 for setLen")
            elseif val < obj.size
                retVal = obj.shorten((obj.size - val));
            elseif val > obj.size
                tFront = obj.front;
                tLen = obj.len;
                tVals = zeros(1, (val - obj.size));
                obj.lengthen(tVals);
                obj.front = tFront;
                obj.len = tLen;
            end
        end

        function retVal = shorten(obj, n)
            arguments
                obj 
                n double = 1 %number of values to shorten 
            end
            newSize = obj.size - n;
            if newSize < 0
                error("shortened queue is smaller than zero")
            end

            toCut = obj.len - newSize;

            retVal = zeros(toCut,1);
            
            for i=1:toCut
                retVal(i) = obj.pop();
            end
            %assign size in rewrite, as the modulo operator in pop() resets the
            %"back" param based on the object size
            obj.rewrite(newSize);
        end

        function val = pop(obj)
            %pop oldest value
            if obj.isempty()
                val = 0;
                return 
            end

            val = obj.vals(obj.back);
            obj.vals(obj.back) = 0;
            obj.len = obj.len - 1;
            obj.back = mod(obj.back, obj.size) + 1;

        end

        function remove(obj, val)
            %remove selected value
            if val == 0
                error("can't remove 0 from the queue")
            end
            obj.rewrite()

            idx = find(obj.vals == val);
            if length(idx) > 1
                error("duplicate values in the queue")
            elseif length(idx) < 1
                warning("the given value was not found in the queue")
                return
            end

            
            obj.vals(idx) = [];
            obj.len = obj.len - 1;
            obj.front = obj.front - 1;
        end
    end
end



%QUEUE CAN DO 4 THINGS
%ADD, assign the next val as what we want, if we reached the queue limit
%then loop back around and overwrite, therefore we first return the current
%value at that idx

%POP, (a specific value, therefore we must reallocate the length of the
%queue)

%LENGTHEN, increase the size of the queue and add new vals at the head

%SHORTEN, shorten the size of the queue and return the last added value

%MAYBE IMPLEMENT IT AS A LINKED LIST INSTEAD OF A QUEUE?