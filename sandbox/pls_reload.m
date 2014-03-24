classdef pls_reload < pls_elem
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        time;
        after_wait;
        ramp_time;
        pos;
    end
    
    methods
        function rd = pls_reload(varargin)
            rd =rd@pls_elem(varargin{:});
        end
        
        function [pulsetab, mktab]=make_tab(rl)
            mktab = zeros(5,0);
            pulsetab = zeros(3, 0);
            if rl.time > 1e-11
                pulsetab(1, 1:4) = cumsum([rl.ramp_time, rl.time, rl.ramp_time,rl.after_wait]);
                pulsetab(2:3, end+(-3:0)) = [repmat(rl.pos', 1, 2), zeros(2)];
            end
        end
    end
    
end

