classdef pls_blank < pls_elem
    %classdef pls_blank < pls_elem
    %   represents a blank pulse element. to be filled in with dictionaries
    %   etc. just requires a name
    
    properties
    end
    
    methods
        function bk = pls_blank(nm)
            bk = bk@pls_elem();
            bk.name = nm;
        end
        
        function [pulsetab, mktab]=make_tab(~)
            pulsetab = zeros(3, 0);
            mktab =  zeros(5, 0);
        end
    end
    
end

