classdef pls_adprep < pls_elem
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function ap = pls_adprep(varargin)
            ap =ap@pls_elem(varargin{:});
            ap.data = struct('time',[],'start',[],'end',[],'mult',1);
        end
        
        function [pulsetab, mktab]=make_tab(ap)
            pulsetab = zeros(3, 0);
            mktab =  zeros(5, 0);
            if ap.data.time > 1e-11
                pulsetab(1, end+(1:2)) = [ap.dt, ap.data.time];
                if length(ap.data.mult)<2
                   dir = ap.data.mult*[1 1]; 
                else
                    dir = ap.data.mult;
                end
                pulsetab(2:3, end-1) = ap.data.start  * dir;
                pulsetab(2:3, end) = ap.data.end  * dir;
            end
            
        end
    end
    
end

