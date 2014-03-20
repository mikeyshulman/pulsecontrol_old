classdef pls_fill < pls_elem
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function rd = pls_fill(varargin)
            rd =rd@pls_elem(varargin{:});
            rd.data = struct('time',[]);
        end
        
        function [pulsetab, mktab]=make_tab(rl)
            mktab = zeros(5,0);
            pulsetab = zeros(3, 0);
        end
    end
    
end

