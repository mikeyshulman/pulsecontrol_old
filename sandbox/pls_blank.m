classdef pls_blank < pls_elem
    %UNTITLED12 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function bk = pls_blank(nm)
            dt = struct('name',nm);
            bk = bk@pls_elem(nm,dt);
        end
        
        function [pulsetab, mktab]=make_tab(~)
            pulsetab = zeros(3, 0);
            mktab =  zeros(5, 0);
        end
    end
    
end

