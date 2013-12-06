classdef sm_pulse < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name;
        data;
        trafofn = '@(x) x';
        pardef = struct('elem_num',{},'par',{},'ind',{});
        param_names={};
        fill_elem = [];
        format = 'elem';
        
    end
    
    methods
        function smp = sm_pulse(stct)
           fn = fieldnames(stct);
           pr = properties('sm_pulse');
           fn = intersect(fn,pr);
           for j = 1:length(fn)
              smp.(fn{j}) = stct.(fn{j}); 
           end
            smp.make_backward_compatible();
           if length(smp.param_names) < length(smp.pardef)
              warning('param names doesn''t have enough entries. this will be confusing...'); 
           end
        end
        
        function make_backward_compatible(smp)
            if isnumeric(smp.pardef)
                pdef = struct();
                for kk = 1:size(smp.pardef,1)
                    pdef(kk).elem_num = smp.pardef(kk,1);
                    if smp.pardef(kk,2)<0
                        pdef(kk).par = 'time';
                    else
                        pdef(kk).par = 'val';
                    end
                    pdef(kk).ind = abs(smp.pardef(kk,2));
                end
                smp.pardef = pdef;
            end
        end
    end
    
end

