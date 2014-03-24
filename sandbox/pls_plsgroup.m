classdef pls_plsgroup < handle
    %sm_pulsegroup < handle: a pulsegroup reprenting an experiment
    %  This is a replacement for a pulsegroup of old plscontrol.
    %  The main difference is that we will treat these plsgrps as
    %  disposible. No logging of group information will occur. All relevant
    %  data for each group will be logged in awgdata and saved with scans.
    %  Groups names could even be left blank if desire (but you shouldn't!)
    %sm_pulsegroup has properies:
        % name: human readable name
        % options: can be 'pack' for speed loading, but not compatible with
            %fancy features like n_reps, jumps
        % pulses: index in plsdata.pulses to take the pulse template from.
            % this can be multivalued to string together different kinds of
            % experiments
        % params: cell array of size plsgrp.pulses to fill the parameters
            % field from the pulse
        % varpar: a cell array, that holds information about the paremeters that
            % will be varied over the group. each element is a struct array
            % with fields:
                % param: which element in params to change
                % val: values that param(param) will take.
                    %needs to be same length for all elements in varpar 
        % trafofn: str or fcn handle to distort wfs, (i.e. rc compensation)
        % chan: which channels of awgdata will get the group
        % readout: struct array with info about readout
        % dict: dictionary to populate pulse elements
        % n_rep: number of times to repeat each pulse (defaults to 1)
        % jump: where each pulse should jump to (defaults to next)
        % pack_data (transient): has wf, marker, tab, data for loading into
            % awg
        % last_update: time of last change to pack_data
    
   % Some userful methods
    % constructor takes a struct with the same fields
    % isequal
    % to_tab: make into table format of just corner points of pulse
    % to_wf: make table into wf
    % make: fill pack_data with all the relevant info
    
                
    properties
        name;
        options = '';
        pulses;
        params;
        varpar;
        trafofn;
        chan;
        readout;
        dict;
        n_rep;
        jump;
    end
    
    properties (SetAccess = immutable)
      type = ''; %human readable flag, for saving/loading
    end
    
    properties (Transient=true)
        last_update;
        %pack_data = struct('tab',{},'marktab',{},'wf',{},'marker',{});
    end
    
    methods
        function pg = pls_plsgroup(config)
            fn = fieldnames(config);
            for j = 1:length(fn)
                if ~strcmp(fn{j},'type')
                    pg.(fn{j}) = config.(fn{j});
                else
                   warning('attempting to set the ''type'' field which can only be set by the constructor. ignoring...'); 
                end
            end
            if isempty(pg.type)
                pg.type = class(pg); % will populate type so we know what class made this object
            end
            pg.make_backward_compatible;
        end
        
        function out = compare(pg,pg2)
            compprops = {'type','pulses','params','varpar','chan','dict'};
            for j = 1:length(compprops)
                out = out && isequal(pg.(compprops{j}),pg2.(compprops{j}));
            end
        end
        
        function b = isequal(pg,pg2)
           p1 = pg.to_struct();
           p2 = pg2.to_struct();
           b = isequal(p1,p2);
        end
        
        function pg=make(plsgrp,ind, clk, tbase, time)
            % ind is the pulses in the group to make
            %pg will be a struct
            %global plsdata;
            if ~exist('time','var')
                time = [];
            end
            if ~exist('tbase','var')|| isempty(tbase)
                global plsdata;
                tbase = plsdata.tbase;
            end
            pg = plsgrp.to_struct();
            if ~iscell(pg.varpar)
                pg.varpar = {pg.varpar};
            end
            if ~iscell(pg.params)
                pg.params = {pg.params};
            end
            if length(pg.pulses)>1
                if length(pg.varpar)==1
                    vp = pg.varpar{1};
                    pg.varpar = cell(1,length(pg.pulses));
                    [pg.varpar]=deal(vp);
                end
                if length(pg.params)==1
                    pp = pg.params{1};
                    pg.params = cell(1,length(pg.pulses));
                    [pg.params] = deal(pp);
                end
            end
            if length(pg.params)~=length(pg.pulses)|| length(pg.varpar)~=length(pg.pulses)
                error('trying to use ''multi'' option with inconsistent dimensions');
            end
            
            % error check the varpars and find out their lengths
            lens = zeros(1,length(pg.varpar));
            if ~isempty(pg.varpar)
                for j = 1:length(pg.varpar)
                    try % hack alert! vertcat is cheap way to check dimensions
                        vertcat(pg.varpar{j}.val);
                        lens(j) = length(pg.varpar{j}(1).val); %the first one will be the right length at this point
                    catch
                        error('varpar element %i has inconsistent dimentions \n',j);
                    end
                end
            end
            
            %lens = cell2mat(cellfun(@(x)length(x),pg.varpar,'UniformOutput',0)); %length of each varpar
            n_pls_tot = max(1,sum(lens));
            %pinds = ones(1,npar); %this will be lookup table to find which varpar and params to use
            pinds = ones(1,lens(1));
            for j = 2:length(lens)
                pinds = [pinds, j*ones(1,lens(j))];
            end
            if ~exist('ind','var') || isempty(ind)
                ind = 1:n_pls_tot;
            end
            
            if isnumeric(pg.pulses)
                global plsdata;
                if any(pg.pulses>length(plsdata.pulses))
                    error('requested pulses %i, but only %i pulses in plsdata',pg.pulses,length(plsdata.pulses))
                end
                pg.pulses = plsdata.pulses(pg.pulses);
            end
            %dp_cp_pls = pg.pulses(pinds(ind));%(ind);
            
            if ~isempty(pg.dict)
                pd=pls_dict.pdpreload(pg.dict,time);
            end
            
            for m = 1:length(ind)
                i = floor((ind(m)-1)/n_pls_tot)+1;
                p = pinds(ind(m));
                j=ind(m)-sum(lens(1:pinds(ind(m))-1));
                
                for jj = 1:length(pg.varpar{p})
                    if ischar(pg.varpar{p}(jj).param)
                       pg.varpar{p}(jj).param = find(ismember(pg.param_names,pg.varpar{p}(jj).param));
                       if isempty(pg.varpar{p}(jj).param)
                           error('no param name found matching param name %i',jj)
                       elseif length(pg.varpar{p}(jj).param)>1
                           error('found more than one param name matching varpar element %i',jj); 
                       end
                    end
                end
                if ~isempty(pg.varpar)
                    vpvals = vertcat(pg.varpar{p}.val); %
                    pparams([pg.varpar{p}.param])=vpvals(:,j);
                    %mask = ~isnan(grpdef.varpar{p}(j, :));
                    %params(end-size(grpdef.varpar{p}, 2) + find(mask)) = grpdef.varpar{p}(j, mask);
                end
                dp_cp_pls = copy(pg.pulses(p));
                if ischar(dp_cp_pls.trafofn)
                    fn = str2func(dp_cp_pls.trafofn);
                else
                   fn = dp_cp_pls.trafofn; 
                end
                pparams = fn(pparams);
                
                % Apply dictionary before varpars; avoids many random bugs.
                if ~isempty(pg.dict)
                    if isstruct(pg.dict)
                        error('not yet implemented');
                        %dp_cp_pls.apply_dict(pg.dict(ind(m)));
                    else
                        dp_cp_pls.apply_dict(pd);
                    end
                    
                end
                
                % update parameters 
                if ~isempty(dp_cp_pls.pardef)
                    pardef = dp_cp_pls.pardef;
                    for n = 1:length(pardef)
                        if isnan(pparams(n))
                            continue;
                        end
                        dp_cp_pls.elems(pardef(n).elem_num).(pardef(n).par)(pardef(n).ind) = pparams(n);
                    end
                end
                pack_data(m) = dp_cp_pls.to_tab;
            end
            %tab_data,clk,tbase,pulse_inds
            pg.pack_data = pls_elem.tabtowf(pack_data,clk,tbase,1:length(ind));
            plsgrp.last_update = now;
        end
        
        function pg_struct = to_struct(pg)
            %this absurd notation with the question mark  pulls class
            %metadata. this will pull the property names and whether they
            %are transient or not (we don't want to save the transient
            %properties. for more help see:
            % http://www.mathworks.com/help/matlab/ref/meta.property.html
            mco = ?pls_plsgroup;
            structify = ~[mco.PropertyList.Transient]; %0 if property is transient
            pp={mco.PropertyList.Name};
            pg_struct = struct();
            for j = 1:length(pp)
                if structify(j)
                    pg_struct.(pp{j}) = pg.(pp{j});
                end
            end
        end
        
        function make_backward_compatible(pg)
            if ~iscell(pg.varpar)
                pg.varpar = {pg.varpar};
            end
            if ~iscell(pg.params)
                pg.params = {pg.params};
            end
            for k = 1:length(pg.varpar)
                if isnumeric(pg.varpar{k}) && ~isempty(pg.varpar{k})
                    vp = struct('param',[],'val',[]);
                    for j = 1:size(pg.varpar{k},2)
                        vp.param(j) = length(pg.params{k})-j+1;
                        vp.val = pg.varpar{k}(:,j);
                    end
                    pg.varpar{k} = vp;
                end
            end
        end
    end %methods
end %class

function [pulse changedout]=pdapply(pd,pulse,time,ind)
%function [pulse changed]=pdapply(pd, pulse,time)
% Apply a pulse dictionary to a pulse.  Return the new pulse.  Changed is
% true if the application was non-trivial
% new feature; an entry like '@1;foo,2;bar,...' will expand to foo for the
%   first pulse, 2 for the second, etc.
if ~exist('time','var')
    time=[];
end
if ~exist('ind','var') %~exist('ind')
    ind='';
end

if 0%~strcmp(pulse.format,'elem')
    changedout=0;
    return;
end

% If pd is a cell array, apply each dictionary in sequence.  This allows
% for some *neat* effects. :p
if iscell(pd)
    changedout=0;
    while 1
        changed = 0;
        for i=1:length(pd)
            [pulse c2] = pdapply(pd{i},pulse,time,ind);
            changed = changed || c2;
        end
        changedout = changed | changedout;
        if ~changed
            break;
        end
    end
    return;
end

if ischar(pd)
    pd=pdload(pd,time);
end
changedout=0;
changed = 1;

while changed
    changed=0;
    for i=1:length(pulse.data)
        if ~changed && strcmp('#',pulse.data(i).type(1))%(pulse.data(i).type(1) == '#')
            entries=regexp(pulse.data(i).type(2:end),'(\d*;)?([^,])*','tokens');
            for e=1:length(entries)
                if isempty(ind) || isempty(entries{e}{1}) || (ind == str2num(entries{e}{1}))
                    changed = 1;
                    pulse.data(i).type = entries{e}{2};
                end
            end
        end
        if ~changed && (pulse.data(i).type(1) == '@')
            if isfield(pd,pulse.data(i).type(2:end))
                %nels=getfield(pd,pulse.data(i).type(2:end));
                nels=pd.(pulse.data(i).type(2:end));
                if ischar(nels)
                    nels={nels};
                end
                if iscell(nels)
                    %nels=struct('type',nels,'time',[],'val',[]);
                    nels=struct('time',[],'type',nels,'val',[]);
                end
                template=pulse.data(i);
                pulse.data(i)=[];
                ot = ~isnan(template.time);
                ov = ~isnan(template.val);
                for j=1:length(nels)
                    if ischar(nels(j))
                        %nels(j) = struct('type',nels(j),'time',[],'val',[]);
                        nels(j) = struct('time',[],'type',nels(j),'val',[]);
                    end
                    nels(j).time(ot)=template.time(ot);
                    nels(j).val(ov)=template.val(ov);
                end
                %pulse.data = [pulse.data(1:i-1) orderfields(nels,pulse.data(1)) pulse.data(i:end)];
                pulse.data = [pulse.data(1:i-1), nels, pulse.data(i:end)];
                %is_hash = [is_hash(1:i),is_hash(i:end)]; %copy the ith is_hash
                %flds = [flds(1:i),flds(i:end)]; %copy the ith flds
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
return
end

% function pardef = bump_pardef(pardef, from, by)
% tobump=find([pardef.elem_num] > from);
% %this for loop is faster and easier to read than the vecotrized version
% for j = tobump
%    pardef(j).elem_num = pardef(j).elem_num+by; 
% end
% %new_nums = num2cell([pardef(tobump).elem_num]+by);
% %[pardef(tobump).elem_num]=deal(new_nums{:});
% end

% function [pd]=pdpreload(pd, time)
% %function [pd]=pdpreload(pd, time)
% % replace any dictionaries in pd with the contents of the dictionary.  Judicous
% % use can vastly speed up pulse dictionaries
% % time, if specified, gives the effective time to load the dictionary at.
% 
% if ~exist('time','var')
%     time=[];
% end
% 
% if iscell(pd)
%    for i=1:length(pd)
%      pd{i} = pdpreload(pd{i},time);
%    end
%    return;
% end
% 
% if ischar(pd)
%     pd=pdload(pd,time);
% end
% end
% 
% function pd = pdload(name, opts)
% %function pd = pdload(name, opts)
% % Load most recent entry in a pulse dictionary.  Load the entire 
% % dictionary if opts='all'
% % if opts is a number, load the most recent dictionary before that.
% %   where opts is a time as returned by 'now' or 'getscantime'
% 
% global plsdata;
% 
% if isstruct(name) % This allows pdload(pload('x')) to be equiv to pdload('x')
%     pd=name;
%     return;
% end
% 
% if ~exist('opts','var')
%     opts = '';
% end
% 
% 
% if isempty(strfind(opts,'all')) && (isempty(opts) || ~isnumeric(opts))
%   load([plsdata.grpdir, 'pd_', name,'_last']);
% else
%   load([plsdata.grpdir, 'pd_', name]);
%   if isnumeric(opts) && ~isempty(opts)
%      times=cellfun(@(x) getfield(x,'time'),pd);
%      i=find(times < opts,1,'last');
%      if isempty(i)
%         error('cannot find dictionary from specified time'); 
%      else
%         pd=pd{i};
%      end
%   end  
% end
%     
% end
