classdef IS_super < matlab.mixin.CustomCompactDisplayProvider
    % Represents an IS with specified ID and ORIENTATION
    
    properties (Constant = true, Hidden = true, Abstract = true)
        POSTFIX  % Reverese, unknown, Forward
    end
    
    properties
        id uint8 = 0;  % IS id
        orientation int8 = 0;  % -1, 0, 1
    end
    
    % Short name (unique)
    properties (Dependent = true)
        name string
    end
    methods
        function name = get.name(obj)
            name = num2str(obj.id) + obj.POSTFIX(obj.orientation + 2);
        end
    end
    
    % Constructor
    methods
        function obj = IS_super(id,orientation)
            if nargin<1
                id = 0;
            end
            if nargin<2
                orientation = 0;
            end
            obj.id = id;
            obj.orientation = orientation;
        end
    end

    % Overloading comparison (EQ, NE) and SORT to allow implementing 
    % UNIQUE and INTERSECT
    methods
        function tf = eq(objs, others)
            objs_name = arrayfun(@(o)o.name, objs);
            others_name = arrayfun(@(o)o.name, others);
            tf = eq(objs_name,others_name);
        end

        function tf = ne(objs, others)
            objs_name = arrayfun(@(o)o.name, objs);
            others_name = arrayfun(@(o)o.name, others);
            tf = ne(objs_name,others_name);
        end

        function [obj_sorted, inds] = sort(objs, varargin)
            [~, inds] = sortrows([[objs.id]', [objs.orientation]'], varargin{:});
            obj_sorted = objs(inds);
        end
    end

    % Compact representation of IS arrays:
    methods (Access = private)
        function str = get_str_representation_for_row_vector(obj)
            if ~isrow(obj)
                error('Only works on rows of ISs')
            end
            strs = arrayfun(@(s) string(sprintf('%3s',s)), [obj.name]);
            str = join(strs,' ');
        end
    end

    methods
        function disp(obj)
            if isrow(obj) && ~isempty(obj)
                str = obj.get_str_representation_for_row_vector();
                disp(str)
            else
                builtin('disp', obj)
            end
        end

        function displayRep = compactRepresentationForSingleLine(obj,displayConfiguration,width)
            displayRep = widthConstrainedDataRepresentation(obj,displayConfiguration,width,...
                StringArray=obj.get_str_representation_for_row_vector());
        end
    end
end
