classdef pulse_dict < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dict; %struct that is a dictionary
        last_update;
        history
    end
    
    methods
        function pd = pulse_dict(s)
           fn = fieldnames(s);
           for j = 1:length(fn)
              pd.(fn{j}) = s.(fn{j}); 
           end
        end
    end
    
end

