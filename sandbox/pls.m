classdef pls
    %classdef pls
    %   terrible way of getting tab completion
    %   put everything in this file that you want to be able to tab
    %   complete into pulse elems. 
    % the to_elem will try to
    %   eval(['pls_',char(obj)]), so things will break if you add to this
    %   list elements that don't exist;
    
    properties
    end
    
    
    enumeration
        reload   %('pls_reload')
        readout  %('pls_readout')
        wait     %('pls_wait')
        adprep   %('pls_adprep')
        adread   %('pls_adread')
        raw      %('pls_raw')
        mark     %('pls_mark')
        start
    end
    
    methods %(Static)
        function e = to_elem(obj)
            e = eval(['pls_',char(obj)]);
        end
    end
    
end

