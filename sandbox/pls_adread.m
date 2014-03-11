classdef pls_adread < pls_elem
    %UNTITLED11 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function ar = pls_adread(varargin)
            ar = ar@pls_elem(varargin{:});
            ar.data = struct('time',[],'st_dly',[],'ed_dly',[],'meas_pt',[],'flag',1);
        end
        
        function [pulsetab, mktab]=make_tab(ar)
            pulsetab = zeros(3, 0);
            mktab =  zeros(5, 0);
            if ar.data.time > 1e-11
                pulsetab(1, end+(1:2)) = [ar.dt, ar.data.time];
                if length(ar.data.mult)<2
                   dir = ar.data.mult*[1 1]; 
                else
                    dir = ar.data.mult;
                end
                pulsetab(2:3, end-1) = ar.data.start  * dir;
                pulsetab(2:3, end) = ar.data.end  * dir;
            end            
        end
    end
    
end

