classdef pls
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
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

