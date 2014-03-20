classdef pls_mark < pls_elem
    %UNTITLED13 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        time;
        val;
    end
    
    methods
        function rw = pls_mark(varargin{:})
            rw.name = 'mark';
        end
        
        function [pulsetab, mktab]=make_tab(mk)
            pulsetab = zeros(3,0);
            mktab =  mk.data.time';
        end
    end
    
end

