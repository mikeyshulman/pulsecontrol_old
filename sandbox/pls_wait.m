classdef pls_wait < pls_elem
    %UNTITLED12 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function wt = pls_wait(varargin)
            wt = wt@pls_elem(varargin{:});
            wt.data = struct('time',[],'val',[],'mult',[]);
        end
        
        function [pulsetab, mktab]=make_tab(wt)
            pulsetab = zeros(3, 0);
            mktab =  zeros(5, 0);
            if wt.data.time > 1e-11
                %fillpos = fillpos +(fillpos==size(pulsetab,2)); % are we filling the wait? if so, we don't want to fill like a ramp
                pulsetab(1, 1:2) = [wt.dt, wt.data.time]; %pinf.tbase*1e6/pinf.clk.
                if ~isempty(wt.data.mult)
                    pulsetab(2:3, end+(-1:0)) = repmat(wt.data.mult*wt.data.val(1:2)', 1, 2);
                else
                    pulsetab(2:3, end+(-1:0)) = repmat(wt.data.val(1:2)', 1, 2);
                end
            end
        end
    end
    
end
