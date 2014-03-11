classdef pls_pls < matlab.mixin.Copyable & handle
    %classdef pls_pls < matlab.mixin.Copyable & handle
    %   represents a pulse template. has properties
    %   elems: array of class pls_elem or derived classes
    %   name:
    %   pardef: struct with fields elem_num, par, ind
    %   param_names
    %   trafofn: string, fn to turn params into pardef
    %   fill: struct with fields elem, time
    %   derives from copyable. allows for pulses to be copied from plsdata
    %   without modifying plsdata;
    
    properties
        elems;
        name;
        trafofn = '@(x) x';
        pardef = struct('elem_num',{},'par',{},'ind',{});
        param_names={};
        fill=struct('elem',{},'time',{});;
    end
    
    methods    
        function pls = pls_pls(stct)
            fn = fieldnames(stct);
            pr = properties('pls_pls');
            fn = intersect(fn,pr);
            for j = 1:length(fn)
                pls.(fn{j}) = stct.(fn{j});
            end
          %  pls.make_backward_compatible();
            if isa(pls.trafofn,'function_handle');
                pls.trafofn = func2str(pls.trafofn);
            end
            if length(pls.param_names) < length(pls.pardef)
                warning('param names doesn''t have enough entries. this will be confusing...');
            end
            for j = 1:length(pls.pardef)
               if isempty(pls.pardef(j).ind) && ~isempty(pls.pardef(j).elem_num)
                  pls.pardef(j).ind = 1; 
               end
            end
        end
        
        function out =to_tab(pls)
            pulsetab = zeros(3, 0);
            mktab =  zeros(5, 0);
            fillpos =[];
            readout = [];
            readpos =[];
            if ~isempty(pls.fill)&&~isempty(pls.fill.elem)
               fill_1 = pls.fill.elem;
            else
                fill_1 = 0;
            end
            fill_2 = [];
            for j = 1:length(pls.elems)
                if isa(pls.elems(j),'pls_fill')
                    fill_2 = [fill_2,j];
                end
            end
            
            if fill_1 && ~isempty(fill_2)
               error('fill_elem given and found pls_fill element'); 
            end
            if length(fill_2)>1
               error('two fills found'); 
            end
            if fill_1
                filltime = pls.fill.time;
                fillelem = pls.fill.elem;
            else
                filltime = pls.elems(fill_2).time;
                fillelem = fill_2;
            end
            
            for j = 1:length(pls.elems)
                if j==fillelem
                   fillpos = size(pulsetab,2);
                   fillmarkpos = size(mktab,2);
                   if isa(pls.elems(j),'pls_wait') || isa(pls.elems(j),'pls_reload')
                      fillpos = fillpos+1; 
                   end
                end
                [t,mt]= pls.elems(j).to_tab;
                if ~isempty(pulsetab)
                    t(1,:) = t(1,:)+pulsetab(1,end);
                end
                if ~isempty(mktab)    
                    mt(1,:) = mt(1,:)+mt(1,end);
                end
                pulsetab = [pulsetab,t];
                mktab = [mktab,mt];
            end
            if ~isempty(fillpos)
                filltime = filltime - pulsetab(1, end);
                if filltime < 0
                    disp(pulsetab);
                    error('Pulse too long by %g (target %g).',-filltime,filltime+pulsetab(1,end));
                end
                pulsetab(1, fillpos+1:end) = pulsetab(1, fillpos+1:end) + filltime;
                if ~isempty(readpos)
                    readout(readpos > fillpos, 2) = readout(readpos > fillpos, 2) + filltime;
                end
                mktab(1, fillmarkpos+1:end) = mktab(1, fillmarkpos+1:end) + filltime;
            end
            out.pulsetab = pulsetab;
            out.marktab = mktab;
            out.readout = readout;
        end
        
        function apply_dict(pulse,pd)
            changed = 1;
            while changed
                changed=0;
                for i=1:length(pulse.elems)
                    if ~changed && isa(pls.elems(i),'pls_blank')
                        if isfield(pd,pulse.elems(i).name)
                            nels=pd.(pulse.elems(i).name);
                            template=pulse.elems(i);
                            ot = ~isnan(template.time);
                            ov = ~isnan(template.val);
                            for j=1:length(nels)
                                fn = fieldname(nels(j).data);
                                for k = 1:length(fn)
                                   if ~isempty(template.data.(fn{k})) && ~isnan(template.data.(fn{k}))
                                      nels(j).data.(fn{k}) = template.data.(fn{k});
                                   end
                                end
                            end
                            pulse.elems = [pulse.elems(1:i-1), nels, pulse.elems(i+1:end)];
                            changed=1;
                            changedout=1;
                            if isfield(pulse,'pardef') && ~isempty(pulse.pardef)
                                pulse.pardef = bump_pardef(pulse.pardef,i,length(nels)-1);
                            end
                            break;
                            
                        end %isfield
                    end
                end
            end
            
        end
    end
    
    methods (Access = protected)
        function cpObj = copyElement(obj)
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            cpObj.elems = copy(obj.elems);
        end
    end
    
    %%%%%
    %%%%%
    methods (Static)
        function new_pardef = make_pardef_compatible(pardef)
            if isstruct(pardef)
                new_pardef = pardef;
                return
            elseif isnumeric(pardef)
                new_pardef = struct();
                for kk = 1:size(pardef,1)
                    new_pardef(kk).elem_num = pardef(kk,1);
                    if pardef(kk,2)<0
                        new_pardef(kk).par = 'time';
                    else
                        new_pardef(kk).par = 'val';
                    end
                    new_pardef(kk).ind = abs(smp.pardef(kk,2));
                end
            else
                error('pardef format not recognized');
            end
        end % end function
    end
    
end

