classdef pls_mark < pls_elem
    %UNTITLED13 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function rw = pls_mark(varargin)
           rw = rw@pls_elem(varargin);
           rw.data.struct('time',[],'val',[]);
           rw.name = mark;
        end
        
        function [pulsetab, mktab]=make_tab(mk)
            pulsetab = zeros(3,0);
            mktab =  mk.data.time';
        end
    end
    
end

