classdef pls_dict < handle & matlab.mixin.Copyable & dynamicprops
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name;
        last_update;
        history;
        time_stamps;
        
    end
    
    methods
        function pd=pls_dict(name,data)
            pd.name = name;
            fn = fieldnames(data);
            for j = 1:length(fn)
               pd.addprop(fn{j});
               pd.(fn{j}) = data.(fn{j});
            end
            
        end
        
        function bool = isequal(pd1,pd2)
           s1 = struct(pd1);
           s2 = struct(pd2);
           ignore_fields = {'history','last_update','time_stamps'};
           s1 = rmfield(s1,ignore_fields);
           s2 = rmfield(s2,ignore_fileds);
           bool = isequal(s1,s2);
        end
        
        function s = to_struct(pd)
            %return all propties that have been added, and also name prop
           pprops = setdiff(properties('pls_dict'),{'name'});
           s = struct(pd);
           s = rmfield(s,pprops);
        end
        
        function update_grps = pdsave(dict,nm)
            global plsdata;
            if exist('nm','var')
                fn = [plsdata.grpdir, 'pd_', nm,'.mat'];
            else
                fn= [plsdata.grpdir, 'pd_', dict.name,'.mat'];
                nm = dict.name;
            end
            if exist(fn,'file')
                try
                    load(fn);            
                    newdict = false;
                catch
                    pd = dict;
                    newdict = true;
                end
            else
                newdict = true;
                pd = dict;
            end
            s = dict.to_struct();
            if isempty(pd.history)
                pd.history=pls_dict(nm,s);
                pd.last_update = now;
                pd.time_stamps = now;
            else
                pd.history(end+1)=pls_dict(nm,dict.to_struct());
                pd.last_update = now;
                pd.time_stamps(end+1) = now;
                fn = fieldnames(s);
                for j = 1:length(fn)
                    if ~isprop(pd,fn{j})
                        pd.addprop(fn{j}); 
                    end
                    pd.(fn{j}) = s.(fn{j});
                end
            end

            if newdict
               fprintf('cannot find file with name %s, saving new dict\n',fn); 
            end
            save(fn,'pd');
            pdlast = copy(pd);
            pdlast.history = [];
            pdlast.time_stamps = [];
            [f1,ext]=strtok(fn,'.'); % break the filename and ".mat"
            save([f1,'_last',ext],'pdlast'); 
            update_grps = [];
            global awgdata
            if ~isempty(awgdata) && ~isempty(awgdata.memory)
               for j = 1:length(awgdata.memory)
                  if any(isequal(awgdata.memory(j).dict) ,nm)
                     update_grps = [update_grps,j];  %#ok<AGROW>
                  end
               end
               if isfield(awgdata,'quiet') && ~awgdata.quiet && changed
                  warning('groups %i use the changed dictionary',update_grps); 
               end
            end
        end
    end
    
    methods (Static)
        function pd = pdload(nm,opts)
            if isa(nm,'pls_dict')
                pd = nm;
                return
            end
            global plsdata;
            if ~exist('opts','var')
                opts = '';
            end
            if isempty(strfind(opts,'all')) && (isempty(opts) || ~isnumeric(opts))
                load([plsdata.grpdir, 'pd_', nm,'_last']);
                pd = pdlast;
            else
                load([plsdata.grpdir, 'pd_', nm]);
                if isnumeric(opts) && ~isempty(opts)
                   ii = find(pd.time_stamps<opts,'last'); 
                   if isempty(ii)
                       error('cannot find dictionary named %s at time %d',pd.name,opts)
                   else
                        pd=pd.history(ii);
                   end
                end
            end    
        end
        
        function [pd]=pdpreload(pd, time)
            %function [pd]=pdpreload(pd, time)
            % replace any dictionaries in pd with the contents of the dictionary.  Judicous
            % use can vastly speed up pulse dictionaries
            % time, if specified, gives the effective time to load the dictionary at.
            
            if ~exist('time','var')
                time=[];
            end
            
            if iscell(pd)
                for i=1:length(pd)
                    pd{i} = pls_dict.pdpreload(pd{i},time);
                end
                return;
            end
            
            if ischar(pd)
                pd=pls_dict.pdload(pd,time);
            end
        end
    end
    
    methods (Access = protected)
        function cpObj = copyElement(obj)
            cpObj = copyElement@matlab.mixin.Copyable(obj);
        end
    end
    
end

