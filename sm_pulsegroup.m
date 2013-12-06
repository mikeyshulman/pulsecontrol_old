classdef sm_pulsegroup < handle
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
                % param_num: which element in params to change
                % val: values that param(param_num) will take.
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
        function pg = sm_pulsegroup(config)
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
        
        function tab_data = totab(pg,pulsedef, p_ind)
            global plsdata;
            %passed the pulsegroup and the index in plsgrp.pulses do turn
            %into a tab. defaults to the first
            % will populate transient variable pack_data
            if nargin < 2
                p_ind = 1;
            end
            dt=-1e-9;  % Shortest meaninful length
            pulsetab = zeros(3, 0);
            mktab =  zeros(5, 0);
            comppos = [];
            fillpos = [];
            rdout = [];
            readpos = [];
            
            % now find the element to fill
            fill_data = find(strcmp(cellstr(char(pulsedef.data.type)),'fill'));
            if isempty(pulsedef.fill_elem) && isempty(fill_data)
                fill_elem = 0; %zero, not empty so can compare to index
            elseif isempty(pulsedef.fill_elem) && ~isempty(fill_data)
                if numel(fill_data)==1
                    fill_elem = fill_data;
                else
                    error('two fill positions found in plsdata.pulses(%i).data \n',plsdata.pulses(p_ind));
                end
            elseif isempty(fill_data) && ~isempty(pulsedef.fill_elem)
                fill_elem = pulsedef.fill_elem;
            else
                error('fill_elem is given even though plsdata.pulses(%i).data has a fill \n',plsdata.pulses(pind));
            end
            
            pulsedef = pulsedef.data;
            
            for i = 1:length(pulsedef)
                if i == fill_elem
                    fillpos = size(pulsetab, 2);
                    filltime = pulsedef(i).time(1);
                    fillmarkpos = size(mktab,2);
                end
                
                switch pulsedef(i).type
                    
                    case 'raw'
                        pulsetab = [pulsetab, [pulsedef(i).time; pulsedef(i).val]];
                        
                    case 'mark'
                        mktab = [mktab, pulsedef(i).time'];
                    case 'fill'
                        fillpos = size(pulsetab, 2);
                        filltime = pulsedef(i).time(1);
                        fillmarkpos = size(mktab,2);      
                    case 'wait'
                        if pulsedef(i).time(1) > 1e-11
                            fillpos = fillpos +(fillpos==size(pulsetab,2)); % are we filling the wait? if so, we don't want to fill like a ramp
                            pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulsedef(i).time(1)]; %pinf.tbase*1e6/pinf.clk.
                            if length(pulsedef(i).val) > 2
                                pulsetab(2:3, end+(-1:0)) = repmat(pulsedef(i).val(3)*pulsedef(i).val(1:2)', 1, 2);
                            else
                                pulsetab(2:3, end+(-1:0)) = repmat(pulsedef(i).val(1:2)', 1, 2);
                            end
                        end
                        
                    case 'reload'
                        % If we're filling the load, push the fillpos 1 forward
                        % so we stretch the wait at the loadpos, not the ramp
                        % to the loadpos
                        % Ignore zero length loads
                        if pulsedef(i).time(2) > 1e-11
                            fillload = (fillpos == size(pulsetab,2));
                            pulsetab(1, end+(1:4)) = pulsetab(1, end) + cumsum(pulsedef(i).time([1 2 1 3]));
                            pulsetab(2:3, end+(-3:0)) = [repmat(pulsedef(i).val(1:2)', 1, 2), zeros(2)];
                            fillpos = fillpos + fillload;
                        end
                    case 'meas_o' % offset measurement
                        pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulsedef(i).time(1)]; %pinf.tbase*1e6/pinf.clk.
                        pulsetab(2:3, end-1) = pulsetab(2:3,end-2);
                        pulsetab(2:3, end) = pulsetab(2:3,end-2);
                        mktab(:, end+1) = [pulsetab(1, end-2)+pulsedef(i).time(2); 0; 0; 0; pulsedef(i).time(1:3)*[1; -1; -1]];
                        if ~isempty(pulsedef(i).val)
                            rdout(end+1, :) = [pulsedef(i).val, pulsetab(1, end-2) + pulsedef(i).time(4), pulsedef(i).time([1 4 5])*[1; -1; -1]];
                            readpos(end+1) = size(pulsetab, 2)-2;
                        end
                        
                    case 'meas'
                        if length(pulsedef(i).val) == 3
                            mpnt = pulsedef(i).val(2:3);
                        else
                            mpnt = [0,0];
                        end
                        pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulsedef(i).time(1)]; %pinf.tbase*1e6/pinf.clk.
                        pulsetab(2:3, end-1) = mpnt;
                        pulsetab(2:3, end) = mpnt;
                        mktab(:, end+1) = [pulsetab(1, end-2)+pulsedef(i).time(2); 0; 0; 0; pulsedef(i).time(1:3)*[1; -1; -1]];
                        if length(pulsedef(i).val) > 0 && ~isnan(pulsedef(i).val(1))
                            rdout(end+1, :) = [pulsedef(i).val(1), pulsetab(1, end-2) + pulsedef(i).time(4), pulsedef(i).time([1 4 5])*[1; -1; -1]];
                            readpos(end+1) = size(pulsetab, 2)-2;
                        end
                        
                    case 'ramp'
                        %allow for multiplies in ramps- helps get direction
                        %right
                        if length(pulsedef(i).val) ==3
                            mult = pulsedef(i).val(3);
                        else
                            mult = 1;
                        end
                        pulsetab(1, end+1) = pulsetab(1, end) + pulsedef(i).time(1);
                        pulsetab(2:3, end) = mult*pulsedef(i).val(1:2);
                        
                    case 'comp'
                        comppos = size(pulsetab, 2)+1;
                        compval  = pulsedef(i).val(1:2);
                        
                        pulsetab(1, end+(1:4)) = pulsetab(1, end) + [0 pulsedef(i).time(2), pulsedef(i).time(1)-sum(pulsedef(i).time(2:3)), ...
                            pulsedef(i).time(1)];
                        pulsetab(2:3, end+(-3:0)) = 0;
                        if length(pulsedef(i).val) >= 4
                            pulsetab(2:3, end) = pulsedef(i).val(3:4);
                        end
                        
                    case 'adprep'
                        if pulsedef(i).time(1) > 1e-11
                            pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulsedef(i).time(1)];
                            if(length(pulsedef(i).val) <= 2)
                                dir=[-1 1];
                            else
                                dir = pulsedef(i).val(3:4);
                            end
                            pulsetab(2:3, end-1) = pulsedef(i).val(1)  * dir;
                            pulsetab(2:3, end) = pulsedef(i).val(2) * dir;
                        end
                    case 'adread'
                        if pulsedef(i).time(1) > 1e-11
                            pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulsedef(i).time(1)];
                            if(length(pulsedef(i).val) <= 2)
                                dir=[-1 1];
                            else
                                dir = pulsedef(i).val(3:4);
                            end
                            pulsetab(2:3, end-1) = pulsedef(i).val(2)  * dir;
                            pulsetab(2:3, end) = pulsedef(i).val(1)  * dir;
                        end
                    case 'RFburst' %added 2013/01/22
                        if length(pulsedef(i).val)==6 && all(~isnan(pulsedef(i).val(4:6)))
                            offst = pulsedef(i).val(4:5);
                            ph = pulsedef(i).val(6);
                        else
                            offst = [0 0];
                            ph = 0;
                        end
                        ts = 0:abs(1e6*dt):pulsedef(i).time(1);
                        if ~exist('firstburst','var') % the first burst will determine the phase for the rest of the pulse
                            firstburst = pulsetab(1,end);
                        end
                        if length(ts) > 1 % number of time points
                            if pulsedef(i).val(1) > 1/abs(dt) % more than 1GHz
                                warning('asking for frequency above 1GHZ');
                            end
                            t_off = pulsetab(1,end)-firstburst;
                            pulsetab(1,end+(1:length(ts)))=pulsetab(1,end)+ts;
                            pulsetab(2,end+(-length(ts)+1:0)) = offst(1)+pulsedef(i).val(2)*sin(ph+2*pi*pulsedef(i).val(1)*(ts+t_off));
                            pulsetab(3,end+(-length(ts)+1:0)) = offst(2)+pulsedef(i).val(3)*sin(ph+2*pi*pulsedef(i).val(1)*(ts+t_off));
                        else
                            warning('asking for only one point of oscillations');
                        end
                    case 'markerburst'
                        ctime = pulsetab(1,end); %current time
                        %first treat it as a wait
                        if pulsedef(i).time(1) > 1e-11
                            fillpos = fillpos +(fillpos==size(pulsetab,2)); % are we filling the wait? if so, we don't want to fill like a ramp
                            pulsetab(1, end+(1:2)) = pulsetab(1, end) + [dt, pulsedef(i).time(1)]; %pinf.tbase*1e6/pinf.clk.
                            if length(pulsedef(i).val) > 2
                                pulsetab(2:3, end+(-1:0)) = repmat(pulsedef(i).val(3)*pulsedef(i).val(1:2)', 1, 2);
                            else
                                pulsetab(2:3, end+(-1:0)) = repmat(pulsedef(i).val(1:2)', 1, 2);
                            end
                            % now deal with the markers
                            %marktab = [marktab, pulsedef(i).time'];
                            if length(pulsedef(i).time)>1 && all(~isnan(pulsedef(i).time(2:3)))
                                stdly = pulsedef(i).time(2); %start delay
                                edly = pulsedef(i).time(3); %end delay
                            else
                                stdly = 0; edly = 0;
                            end
                            mk = pulsedef(i).val(4:7); mk(isnan(mk))=0;
                            mktab = [mktab, [ctime-stdly, (pulsedef(i).time(1)+edly+stdly)*mk]'];
                        end
                    otherwise
                        error('Invalid pulse element %i: %s.\n', i, pulsedef(i).type)
                end
            end
            
            if ~isempty(comppos)
                pulsetab(2:3, comppos+(1:2)) = 2*repmat(compval(1:2)', 1, 2);
            end
            
            %pulsetab(2:3, :) = pulsetab(2:3, :)./pinf.scale;
            %pinf = rmfield(pinf, 'scale');
            
            if ~isempty(fillpos)
                filltime = filltime - pulsetab(1, end);
                if filltime < 0
                    pulsetab
                    error('Pulse too long by %g (target %g).',-filltime,filltime+pulsetab(1,end));
                end
                pulsetab(1, fillpos+1:end) = pulsetab(1, fillpos+1:end) + filltime;
                if ~isempty(readpos)
                    rdout(readpos > fillpos, 2) = rdout(readpos > fillpos, 2) + filltime;
                end
                mktab(1, fillmarkpos+1:end) = mktab(1, fillmarkpos+1:end) + filltime;
            end
            
            mask = all(abs(diff(pulsetab(2:3, :), [], 2)) < 1e-14);
            pulsetab(:, [false, mask(2:end)&mask(1:end-1)]) = [];
            tab_data.marktab = sortrows(mktab',1);
            tab_data.readout = rdout;
            tab_data.tab = pulsetab;
        end
        
        function tab_data=towf(pg,tab_data,clk,tbase,pulse_inds)
            dt=1e-11;
            %looks at pack_data(pulse_inds).tab and marktab and makes wf
            for jj = pulse_inds
                if isempty(tab_data(jj).tab)||isempty(tab_data(jj).marktab)
                    error('pulse number %i not filled into table\n',jj)
                end
                if 0
                    %             % these fields are optional, no need to store in database
                    %             if ~isfield(pulseinf, 'marktab')
                    %                 pulseinf.marktab = [];
                    %             end
                    %
                    %             if ~isfield(pulseinf, 'pulsefn')
                    %                 pulseinf.pulsefn = [];
                    %             elseif ~isempty(pulseinf.pulsefn) && ~isfield(pulseinf.pulsefn, 'args')
                    %                 [pulseinf.pulsefn.args] = deal(cell(2, 0));
                    %             end
                    %
                    %             if ~isfield(pulseinf, 'readout')
                    %                 pulseinf.readout = [];
                    %             end
                end
                pulsetab=tab_data.tab;
                nchan = size(pulsetab, 1)-1;
                npoints = round(max(pulsetab(1, :)) * tbase * clk/1e9);
                
                data = zeros(nchan, npoints+1);
                ttime = linspace(pulsetab(1, 1), pulsetab(1, end), npoints+1);
                
                %            avg = zeros(nchan, 1);
                for j = 1:nchan
                    for i = 2:size(pulsetab, 2)
                        mask = ttime >= pulsetab(1, i-1)-dt & ttime <= pulsetab(1, i)+dt;
                        % added small shifts to mitigate rounding errors 08/04/09. Never seen to matter.
                        % below makes writes the pulse into data using lines to connect the
                        % corners defined in pulstab
                        if 0
                            data(j, mask) = (-pulsetab(j+1, i-1) * (ttime(mask) - pulsetab(1, i)) ...
                                + pulsetab(j+1, i) * (ttime(mask) - pulsetab(1, i-1)))./...
                                (pulsetab(1, i) -  pulsetab(1, i-1));
                        else
                            data(j, mask) = ((-pulsetab(j+1, i-1) + pulsetab(j+1,i)) * ttime(mask) + ...
                                pulsetab(j+1,i-1) * pulsetab(1, i) - pulsetab(j+1,i) * pulsetab(1, i-1))./...
                                (pulsetab(1, i) -  pulsetab(1, i-1));
                        end
                        
                    end
                end
                % lets pulses be defined with functions (eg. sin, cos) instead of
                % just lines
                %             for i = 1:length(pulseinf.pulsefn)
                %                 mask = ttime > pulseinf.pulsefn(i).t(1) & ttime <= pulseinf.pulsefn(i).t(2);
                %                 for j = 1:nchan
                %                     data(j, mask) = pulseinf.pulsefn(i).fn{j}(ttime(mask)-pulseinf.pulsefn(i).t(1), pulseinf.pulsefn(i).args{j, :}) - avg(j);
                %                 end
                %             end
                
                data(:, end) = [];
                marker = zeros(nchan, npoints, 'uint8');
                
                % extend marktab to be right dimensions
                mktab = tab_data(jj).marktab;
                mktab(end+1:2*nchan+1,:) = 0;
                
                for i = 1:size(mktab, 2)
                    for j = 1:nchan
                        for k = 1:2;
                            mask = ttime(1:end-1) >= mktab(1, i) - dt &...
                                ttime(1:end-1) < mktab(1, i) + mktab(j*2+k-1, i)-2e-11;
                            marker(j, mask) = bitor(marker(j, mask), k);
                        end
                    end
                end
                %plsgrp.pack_data(jj).marktab = mktab;
                tab_data(jj).marker = marker;
                tab_data(jj).wf = data;%+vc;
                %plsgrp.readout = pulseinf.readout;
                plsgrp.last_update=now;
            end
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
            
            pg.pulses = plsgrp.default(pg.pulses);
            plsdef = pg.pulses(pinds(ind));%(ind);
            
            if ~isempty(pg.dict)
                pg.dict=pdpreload(pg.dict,time);
            end
            pg=rmfield(pg,'pulses');
            
            for m = 1:length(ind)
                i = floor((ind(m)-1)/n_pls_tot)+1;
                p = pinds(ind(m));
                j=ind(m)-sum(lens(1:pinds(ind(m))-1));
                
                pparams = pg.params{p};
                if ~isempty(pg.varpar)
                    vpvals = vertcat(pg.varpar{p}.val); %
                    pparams([pg.varpar{p}.param_num])=vpvals(:,j);
                    %mask = ~isnan(grpdef.varpar{p}(j, :));
                    %params(end-size(grpdef.varpar{p}, 2) + find(mask)) = grpdef.varpar{p}(j, mask);
                end
                
                if ischar(plsdef(ind(m)).trafofn)
                    fn = str2func(plsdef(ind(m)).trafofn);
                else
                   fn = plsdef(ind(m)).trafofn; 
                end
                pparams = fn(pparams);
                
                % Apply dictionary before varpars; avoids many random bugs.
                if ~isempty(pg.dict)
                    if isstruct(pg.dict)
                        plsdef(m) = pdapply(pg.dict(ind(m)),plsdef(m),time,ind(m));
                    else
                        plsdef(m) = pdapply(pg.dict,plsdef(m),time,ind(m));
                    end
                    
                end
                
                %                mask = ~isnan(pparams);
                % update parameters - could move to plstowf
                if ~isempty(plsdef(m).pardef)
                    pardef = plsdef(m).pardef;
                    for n = 1:length(pardef)
                        if isnan(pparams(n))
                            continue;
                        end
                        plsdef(m).data(pardef(n).elem_num).(pardef(n).par)(pardef(n).ind) = pparams(n);
                    end
                end
                pack_data(m) = plsgrp.totab(plsdef(m),ind(m));
                % towf(pg,clk,tbase,pulse_inds)
                %plsgrp.towf(tab_data,clk,tbase, ind(m));
                %plsgrp.last_update = now;
            end
            pg.pack_data=plsgrp.towf(pack_data,clk,tbase,1:length(ind));
        end
        
        function pulse = default(~,plsnum)
            global plsdata;
            if plsnum > length(plsdata.pulses)
                error('Requested pulse %d, but only %d are defined.  Did you plssync?',plsnum,length(plsdata.pulses));
            end
            pulse = plsdata.pulses(plsnum);
            return;
        end
        
        function pg_struct = to_struct(pg)
            %this absurd notation with the question mark  pulls class
            %metadata. this will pull the property names and whether they
            %are transient or not (we don't want to save the transient
            %properties. for more help see:
            % http://www.mathworks.com/help/matlab/ref/meta.property.html
            mco = ?sm_pulsegroup;
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
                if isnumeric(pg.varpar{k}) && ~ismepty(pg.varpar{k})
                    vp = struct('param_num',{},'val',{});
                    for j = 1:size(pg.varpar{k},2)
                        vp.param_num(j) = length(pg.params{k})-j+1;
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
        if ~changed && (pulse.data(i).type(1) == '#')
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
                    nels=struct('type',nels,'time',[],'val',[]);
                end
                template=pulse.data(i);
                pulse.data(i)=[];
                ot = ~isnan(template.time);
                ov = ~isnan(template.val);
                for j=1:length(nels)
                    if ischar(nels(j))
                        nels(j) = struct('type',nels(j),'time',[],'val',[]);
                    end
                    nels(j).time(ot)=template.time(ot);
                    nels(j).val(ov)=template.val(ov);
                end
                pulse.data = [pulse.data(1:i-1) orderfields(nels,pulse.data(1)) pulse.data(i:end)];
                changed=1;
                changedout=1;
                if isfield(pulse,'pardef') && ~isempty(pulse.pardef)
                    pulse.pardef = bump_pardef(pulse.pardef,i,length(nels)-1);
                end
                break;
                
            end
        end
    end
end
return
end

function pardef = bump_pardef(pardef, from, by)
tobump=find([pardef.elem_num] > from);
%FIXME this is stupid
% for j = tobump
%    pardef(j).elem_num = pardef(j).elem_num+by; 
% end
new_nums = num2cell([pardef(tobump).elem_num]+by);
[pardef(tobump).elem_num]=deal(new_nums{:});
%[pardef(tobump).elem_num]=deal(pardef(tobump).elem_num);
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
     pd{i} = pdpreload(pd{i},time);
   end
   return;
end

if ischar(pd)
    pd=pdload(pd,time);
end
end

function pd = pdload(name, opts)
%function pd = pdload(name, opts)
% Load most recent entry in a pulse dictionary.  Load the entire 
% dictionary if opts='all'
% if opts is a number, load the most recent dictionary before that.
%   where opts is a time as returned by 'now' or 'getscantime'

% (c) 2010 Hendrik Bluhm.  Please see LICENSE and COPYRIGHT information in plssetup.m.

global plsdata;

if isstruct(name) % This allows pdload(pload('x')) to be equiv to pdload('x')
    pd=name;
    return;
end

if ~exist('opts','var')
    opts = '';
end


if isempty(strfind(opts,'all')) && (isempty(opts) || ~isnumeric(opts))
  load([plsdata.grpdir, 'pd_', name,'_last']);
else
  load([plsdata.grpdir, 'pd_', name]);
  if isnumeric(opts) && ~isempty(opts)
     times=cellfun(@(x) getfield(x,'time'),pd);
     i=find(times < opts,1,'last');
     if isempty(i)
        error('cannot find dictionary from specified time'); 
     else
        pd=pd{i};
     end
  end  
end
    
end
