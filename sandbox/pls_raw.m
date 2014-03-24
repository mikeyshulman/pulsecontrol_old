classdef pls_raw < pls_elem
    %UNTITLED13 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        time;
        val;
    end
    
    methods
        function rw = pls_raw(varargin)
           rw = rw@pls_elem(varargin);
           rw.name = 'raw';
        end
        
        function [pulsetab, mktab]=make_tab(rw)
            pulsetab = [rw.time;rw.val];
            mktab =  zeros(5, 0);
        end
    end
    
end

