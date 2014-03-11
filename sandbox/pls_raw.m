classdef pls_raw < pls_elem
    %UNTITLED13 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function rw = pls_raw(varargin)
           rw = rw@pls_elem(varargin);
           rw.data=struct('time',[],'val',[]);
           rw.name = 'raw';
        end
        
        function [pulsetab, mktab]=make_tab(rw)
            pulsetab = [rw.data.time;rw.data.val];
            mktab =  zeros(5, 0);
        end
    end
    
end

