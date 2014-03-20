classdef pls_readout < pls_elem
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function rd = pls_readout(varargin)
            rd =rd@pls_elem(varargin{:});
            rd.data = struct('time',1,'st_dly',.2,'ed_dly',0,'meas_pt',[],'flag',1);
        end
        
        function [pulsetab, mktab] = make_tab(rd)
            pulsetab = zeros(3, 0);
            mktab =  zeros(5, 0);
            if isempty(rd.data.meas_pt)
                rd.data.meas_pt = [0 0];
            end
            pulsetab(1, end+(1:2)) = [rd.dt, rd.data.time]; %pinf.tbase*1e6/pinf.clk.
            pulsetab(2:3, end-1) = rd.data.meas_pt;
            pulsetab(2:3, end) = rd.data.meas_pt;
            mktab(:, 1) = [rd.data.st_dly; 0; 0; 0; [rd.data.time-rd.data.st_dly-rd.data.ed_dly]];
        end
    end
    
end

