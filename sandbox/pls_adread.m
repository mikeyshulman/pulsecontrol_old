classdef pls_adread < pls_elem
    %UNTITLED11 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        time;
        start;
        finish;
        mult;
    end
    
    methods
        function ar = pls_adread(varargin)
            ar = ar@pls_elem(varargin{:});
        end
        
        function [pulsetab, mktab]=make_tab(ar)
            pulsetab = zeros(3, 0);
            mktab =  zeros(5, 0);
            if ar.time > 1e-11
                pulsetab(1, end+(1:2)) = [ar.dt, ar.time];
                if length(ar.mult)<2
                   dir = ar.mult*[1 1]; 
                else
                    dir = ar.mult;
                end
                pulsetab(2:3, end-1) = ar.start  * dir;
                pulsetab(2:3, end) = ar.finish  * dir;
            end            
        end
    end
    
end

