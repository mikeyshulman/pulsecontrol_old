classdef sm_awg < handle
    %sm_awg < handle: represents an awg
    %   sm_awg should usually be accessed through sm_awg_data. 
    %   the sm_awg class has properties:
            % name: mostly useless
            % inst: tcpip, visa, gbip, etc. for communication
            % chans: channel numbers, i.e. 1,2,3,4
            % scale: the Vpp of each channel
            % clk: the clk frequency:
            % offst: offset in mV for each channel
            % bits: voltage resolution
            % slave: 0 if master. for triggering one off another
            % triglen: length of trigger pls
            % trig_pls_name:
            % zero_pls: struct with fields
                % zero_chan: index for each channel of which zeropls to
                    % use. this will be different for channels with different
                    % offsets
                % zero_lens: cell array size unique(zero_chan) x N where N
                    % is the number of different length pulses are used
            %type: i.e. '5k'
            % pulsegroups: struct array of pulsegroups with relevant info
                % for loading them and running scans with them
            % wf_dir: directory where wfs and .seq files go
            % waveforms: a length(chans) x n_groups cell array with the
                % names of waveforms for loading
      
    % methods should be called from awg_data to avoid confounding the user
    % constructor takes a struct with the same fields
    
    properties
        inst;
    end
    properties (SetAccess = {?sm_awg, ?sm_awg_data})
        name = 'AWG5014C';
        chans = [1 2 3 4];
        scale = [1 1 1 1]*600; %Vpp
        clk = 1e9;
        is5k =1; %is it an awg5k (for raw/amp etc)
        offst = [0 0 0 0];
        bits = 14;
        slave = 0;
        trig_pls = struct('len',1000,'name','trig_awg1','chan',2,'marker',1);
        zero_pls = struct('zero_chan',[],'zero_name','zero_%d_%04d'); % initialize to empty so constructor can look and populate
        pulsegroups;
        wf_dir;
        waveforms = {}; % names of the wfs
        pls_lens= []; %keep track of pulse length
    end
    properties (SetAccess = immutable)
       type=''; 
    end
    
    methods
        function awg = sm_awg(s)
            fn = fieldnames(s);
            for j = 1:length(fn)
                if ~strcmp(fn{j},'type')
                    awg.(fn{j}) = s.(fn{j});
                else
                    warning('attempting to set the ''type'' field which can only be set by the constructor. ignoring...');
                end
            end
            if isempty(awg.type)
                awg.type = class(awg); % will populate type so we know what class made this object
            end
            if 1%isempty(awg.zero_pls)
                [~,~,inds]= unique(awg.offst);
                awg.zero_pls(1).zero_chan = inds;
            end
        end
        
         function rm_seq(awg,grp,ctrl) % remove groups from seq table
%             if strcmp(grp, 'all')
%                 awg.stop();
%                 fprintf(awg.inst,'SEQ:LENG 0');
%                 awg.pulsegroups= {};
%                 awg.seqpulses = [];
%                 return;
%             end
%             
%             if ~isnumeric(grp)
%                 grp=strmatch(grp, strvcat(awgdata.pulsegroups.name), 'exact');
%                 %grp = awggrpind(grp); 
%             end
%             grp(2:end) = [];
%             if isnan(grp)
%                 return;
%             end
%             awg.stop();
%             
%             if exist('ctrl','var') && strfind(ctrl, 'after')
%                 fprintf(awg.inst, 'SEQ:LENG %d', awg.pulsegroups(grp).seqind-1 + sum(awgdata(a).pulsegroups(grp).nline));
%                 awg.seqpulses(awg.pulsegroups(grp).seqind + sum(awgdata(a).pulsegroups(grp).npulse):end) = [];
%                 awg.pulsegroups(grp+1:end) = [];
%                 % may miss trigger line.
%                 return;
%             end
%             
%             fprintf(awg.inst, 'SEQ:LENG %d', awg.pulsegroups(grp).seqind-1);
%             awg.seqpulses(awgdata(a).pulsegroups(grp).seqind:end) = [];
%             groups = {awg.pulsegroups(grp+1:end).name};
%             awg.pulsegroups(grp:end) = [];     
         end
        
        function add_trig_wf(awg)
           return; %FIMXE 
        end
        
        function lint(awg)
           if isempty(awg.wf_dir)
               warning('wf_dir is empty. ain''t gonna be easy to write waveforms\n');
           end
        end
                
         function add(awg,groups)
%             dosave = false; % keeps track of whether awgdata changed.
%             gstart=toc;
%             for k = 1:length(groups)
%                 plsgrp = sm_loadplsgrp([plsdata.grpdir, 'pg_', groups{k}]);
%                 
%                 while plsgrp.isstale()
%                     if ~awg(1).quiet
%                         fprintf('Latest pulses of group %s not loaded; %s > %s.\n', groups{k}, ...
%                             datestr(lastupdate), datestr(plslog(end).time(end)));
%                     end
%                     tstart=toc;
%                     plsmakegrp(groups{k},'upload');
%                     ts2 = toc;
%                     awgcntrl('wait');
%                     %  fprintf('Load time=%f secs, wait time=%f\n',toc-tstart,toc-ts2);
%                     load([plsdata.grpdir, 'pg_', groups{k} '.mat']);
%                 end
%                 
%                 
%                 if strcmp(grpdef.ctrl(1:min([end find(grpdef.ctrl == ' ', 1)-1])), 'grp')...
%                         && ~isempty(strfind(grpdef.ctrl, 'seq')) % combine groups at sequence level.
%                     
%                     % retrieve channels of component groups
%                     clear chan;
%                     lastload = 0;
%                     for m = 1:length(grpdef.pulses.groups)
%                         gd=plsinfo('gd', grpdef.pulses.groups{m});
%                         pl=plsinfo('log', grpdef.pulses.groups{m});
%                         if pl(end).time > lastload
%                             lastload = pl(end).time;
%                         end
%                         rf={'varpar'}; % Required fields that may be missing
%                         for qq=1:length(rf) % Put in empty copies of all required fields.
%                             if ~isfield(gd,rf{qq})
%                                 gd=setfield(gd,rf{qq},[]);
%                             end
%                         end
%                         chan(m) = orderfields(gd);
%                     end
%                     chan = {chan.chan};
%                     seqmerge = true;
%                 else
%                     if ~isfield(grpdef, 'pulseind')
%                         grpdef.pulseind = 1:size(zerolen, 1);
%                     end
%                     
%                     seqmerge = false;
%                 end
%                 
%                 if ~isfield(grpdef, 'nrep')
%                     grpdef.nrep = 1;
%                 end
%                 
%                 for a=1:length(awgdata)
%                     npls = size(grpdef.pulseind, 2);
%                     nchan = length(awgdata(a).chans); % alternatively use awgdata or data size
%                     usetrig = (grpdef.nrep(1) ~= Inf) && isempty(strfind(grpdef.ctrl, 'notrig'));
%                     
%                     if isempty(awgdata(a).pulsegroups)
%                         startline = 1;
%                         gind = [];
%                     else
%                         gind = strmatch(grpdef.name, {awgdata(a).pulsegroups.name}, 'exact');
%                         if ~isempty(gind)
%                             startline = awgdata(a).pulsegroups(gind(1)).seqind;
%                             if npls + usetrig ~= sum(awgdata(a).pulsegroups(gind(1)).npulse);
%                                 error('Number of pulses changed in group %s. Use awgrm first!', grpdef.name);
%                             end
%                         else
%                             startline = awgdata(a).pulsegroups(end).seqind + sum(awgdata(a).pulsegroups(end).nline);
%                         end
%                     end
%                     
%                     if ~exist('zerolen','var')  % For sequence combined groups, getting this is hard.
%                         zerolen = plsinfo('zl',grpdef.name);
%                     end
%                     if isempty(gind) % group not loaded yet, extend sequence
%                         
%                         gind = length(awgdata(a).pulsegroups)+1;
%                         awgdata(a).pulsegroups(gind).name = grpdef.name;
%                         awgdata(a).pulsegroups(gind).seqind = startline;
%                         awgdata(a).pulsegroups(gind).lastupdate = lastupdate;
%                         awgdata(a).pulsegroups(gind).npulse = [npls usetrig];
%                         if strfind(grpdef.ctrl,'pack')
%                             awgdata(a).pulsegroups(gind).nline = 1+usetrig;
%                             % Hack alert; way too much code assumes zl == pulselen.  For
%                             % packed groups, we ignore it and work out the correct length
%                             % ourselves.
%                             zlmult=npls;
%                             npls=1;
%                         else
%                             zlmult=1;
%                             awgdata(a).pulsegroups(gind).nline = npls+usetrig;
%                         end
%                         fprintf(awgdata(a).awg, sprintf('SEQ:LENG %d', startline + awgdata(a).pulsegroups(gind).nline-1));
%                         dosave = 1;
%                     else
%                         if strfind(grpdef.ctrl,'pack')
%                             zlmult = npls;
%                             npls=1;
%                         else
%                             zlmult=1;
%                         end
%                         if any(awgdata(a).pulsegroups(gind).nrep ~= grpdef.nrep) || ...
%                                 (isfield(plslog(end),'readout') && (any(any(awgdata(a).pulsegroups(gind).readout ~= plslog(end).readout)))) || ...
%                                 any(any(awgdata(a).pulsegroups(gind).zerolen ~= zerolen)) % nrep or similar changed
%                             dosave = 1;
%                         end
%                     end
%                     
%                     if ~isfield(grpdef, 'jump')
%                         if strfind(grpdef.ctrl, 'loop')
%                             grpdef.jump = [npls; 1];
%                         else
%                             grpdef.jump = [];
%                         end
%                     end
%                     
%                     
%                     % Cache frequently used information in awgdata
%                     awgdata(a).pulsegroups(gind).nrep = grpdef.nrep;
%                     if seqmerge
%                         awgdata(a).pulsegroups(gind).lastload = lastload;
%                     else
%                         awgdata(a).pulsegroups(gind).lastload = plslog(end).time(1);
%                     end
%                     
%                     awgdata(a).pulsegroups(gind).zerolen = zerolen;
%                     if isfield(plslog(end),'readout')
%                         awgdata(a).pulsegroups(gind).readout=plslog(end).readout;
%                     else
%                         awgdata(a).pulsegroups(gind).readout = plsinfo('ro',grpdef.name);
%                     end
%                     
%                     if usetrig
%                         for j = 1:nchan
%                             if mod(j,2) == 1
%                                 fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "trig_%08d"', startline, j, awgdata(a).triglen));
%                             else
%                                 fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "zero_%08d_%d"', startline, j, awgdata(a).triglen, awgdata(a).zerochan(j)));
%                             end
%                         end
%                         if isfield(awgdata(a),'slave') && ~isempty(awgdata(a).slave) && (awgdata(a).slave)
%                             fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:TWAIT 1\n', startline));
%                         end
%                     end
%                     
%                     
%                     for i = 1:npls
%                         ind = i-1 + startline + usetrig;
%                         if ~seqmerge % pulses combined here.
%                             for j = 1:nchan
%                                 ch = find(awgdata(a).chans(j) == grpdef.chan);
%                                 if ~isempty(ch) &&  zerolen(grpdef.pulseind(i), ch) < 0
%                                     % channel in group and not zero
%                                     fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "%s_%05d_%d"', ind, j, ...
%                                         grpdef.name, grpdef.pulseind(i), ch));
%                                 else
%                                     % hack alert. We should really make zerolen a cell array.  fixme.
%                                     fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "zero_%08d_%d"', ind, j, ...
%                                         zlmult*abs(zerolen(grpdef.pulseind(i), 1))*awgdata(a).clk/awgdata(1).clk,awgdata(a).zerochan(j)));
%                                 end
%                             end
%                         else
%                             for m = 1:length(grpdef.pulses.groups)
%                                 for j = 1:length(chan{m}) % channels of component groups
%                                     ch = find(awgdata(a).chans == chan{m}(j));
%                                     if ~isempty(ch)
%                                         %if 1 % zero replacement not implemented
%                                         fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "%s_%05d_%d"', ind, ch, ...
%                                             grpdef.pulses.groups{m}, grpdef.pulseind(m, i), ch));
%                                     else
%                                         error('This won''t work');
%                                     end
%                                     %else
%                                     %fprintf(awgdata.awg, sprintf('SEQ:ELEM%d:WAV%d "zero_%08d"', ind, awgdata.chans(j), ...
%                                     %    abs(zerolen(grpdef.pulseind(i), 1))));
%                                     %end
%                                 end
%                             end
%                         end
%                         if grpdef.nrep(min(i, end)) == Inf  || grpdef.nrep(min(i, end)) == 0 ...
%                                 || (i == npls && isempty(strfind(grpdef.ctrl, 'loop')) && (isempty(grpdef.jump) || all(grpdef.jump(1, :) ~= i)))
%                             fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:LOOP:INF 1', ind));
%                         else
%                             fprintf(awgdata(a).awg, 'SEQ:ELEM%d:LOOP:INF 0', ind); % default
%                             fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:LOOP:COUN %d', ind, grpdef.nrep(min(i, end))));
%                         end
%                         
%                         fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:GOTO:STAT 0', ind));
%                         
%                         if grpdef.nrep(min(i, end)) == Inf && isreal(grpdef.pulses) &&  ...
%                                 (length(awgdata(a).seqpulses) < ind || awgdata(a).seqpulses(ind) ~= grpdef.pulses(grpdef.pulseind(i)));
%                             dosave = 1;
%                             awgdata(a).seqpulses(ind) = grpdef.pulses(grpdef.pulseind(i));
%                         end
%                         if ~mod(i, 100) && ~awgdata(1).quiet
%                             fprintf('%i/%i pulses added.\n', i, npls);
%                         end
%                     end
%                     %fprintf('Group load time: %g secs\n',toc-gstart);
%                     
%                     jstart=toc;
%                     % event jumps
%                     %SEQ:ELEM%d:JTARget:IND
%                     %SEQ:ELEM%d:JTARget:TYPE
%                     
%                     for j = 1:size(grpdef.jump, 2)
%                         fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:GOTO:IND %d', startline+usetrig-1 + grpdef.jump(:, j)));
%                         fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:GOTO:STAT 1', startline+usetrig-1 + grpdef.jump(1, j)));
%                     end
%                 end
%                 if ~exist('seqlog','var')
%                     seqlog.time = now;
%                 else
%                     seqlog(end+1).time = now;
%                 end
%                 seqlog(end).nrep = grpdef.nrep;
%                 seqlog(end).jump = grpdef.jump;
%                 
%                 save([plsdata.grpdir, 'pg_', groups{k}], '-append', 'seqlog');
%                 %fprintf('Jump program time: %f secs\n',toc-jstart);
%                 wstart=toc;
%                 awgcntrl('wait');
%                 %fprintf('Wait time: %f secs; total time %f secs\n',toc-wstart,toc-astart);
%                 nerr=0;
%                 for a=1:length(awgdata(a))
%                     err=query(awgdata(a).awg, 'SYST:ERR?');
%                     if ~isempty(strfind(err, 'No error'))
%                         nerr=nerr+1;
%                     end
%                 end
%                 if nerr == 0
%                     if ~awgdata(1).quiet
%                         fprintf('Added group %s on index %i. %s', grpdef.name, gind, err);
%                     end
%                     logentry('Added group %s on index %i.', grpdef.name, gind);
%                 end
%             end
%             if dosave
%                 awgsavedata;
%             end
%             
%             
         end
        
        function out = is_wf_list_empty(awg)
          out= query(awg.inst, 'WLIS:SIZE?', '%s\n', '%i') == 25; % nothing loaded (except predefined) 
        end
        
        function was_synced = syncwaveforms(awg)
            % Make sure the list of pulses is awgdata is consistent with the awg.
            % we assume if the number of pulses is right, everything is.\
            awg.clr();
            npls=str2num(query(awg.inst,'WLIS:SIZE?'));
            if isfield(awgdata(a),'waveforms') && (length(awg.waveforms) == npls)
                was_synced = 1;
                return;
            end
            fprintf('AWG waveform list out of date.  Syncing.');
            awg.waveforms=cell(npls,1);
            for l=1:npls
                r=query(awg.inst,sprintf('WLIS:NAME? %d',l-1));
                awgdata(a).waveforms{l}=r(2:end-2);
            end
            fprintf('..  Done.\n');
            was_synced = 0;
        end
        
        function zerolen = loadgrp(grp, ind)
            
            awg.syncwaveforms(); % make sure the waveform list is up-to-date.
            dosave = false;
            dind=find([grp.pulses(1).data.clk] == awg.clk); 
            nchan = length(awg.chans); % alternatively use awgdata or data size

            % fixme; emit an error if this changes zerochan and wlist.size ~= 25.
            % The screwiness here is to get each channel with a unique offset/scale combo
            [offsets offsetchan awgdata(a).zerochan] = unique(awg.offset./awg.scale);
            offsets=awg.offset(offsetchan);
            
            % create trig pulse (and corresponding 0) if waveform list empty.
            if awg.is_wf_list_empty() % nothing loaded (except predefined)
                zdata=zeros(1,awgdata(a).triglen);
                zmarker=repmat(1,1,awgdata(a).triglen); % channel 1&3 marker 1
                awgloadwfm(a,zdata,zmarker,sprintf('trig_%08d',awg.triglen),1,1);
                zmarker=zeros(1,awgdata(a).triglen);
                for l=1:length(offsets)
                    awg.load_raw_wf(a,zdata,zmarker,sprintf('zero_%08d_%d',awg.triglen,l),offsetchan(l),1);
                end
                dosave = 1;
                awg.zero_pls_lengths = awg.triglen;
            end
            nzpls=0;
            for i = 1:length(grp.pulses)
                npts = size(grp.pulses(i).data(dind).wf, 2);
                if ~any(awgdata(a).zeropls == npts) % create zero if not existing yet
                    zdata=zeros(1,npts);
                    for l=1:length(offsets)
                        zname=sprintf('zero_%08d_%d', npts, l);
                        awgloadwfm(a,zdata,zdata,zname,offsetchan(l),1);
                    end
                    zdata=[];
                    awg.zero_pls_lengths(end+1) = npts;
                    dosave = 1;
                end
                
                for j = 1:size(grp.pulses(i).data(dind).wf, 1)
                    ch=find(grp.chan(j)==awgdata(a).chans);
                    if isempty(ch)
                        continue;
                    end
                    if any(abs(grp.pulses(i).data(dind).wf(j, :)) > awgdata(a).scale(ch)/(2^awgdata(a).bits)) || any(grp.pulses(i).data(dind).marker(j,:) ~= 0)
                        name = sprintf('%s_%05d_%d', grp.name, ind(i), j);
                        
                        if isempty(strmatch(name,awgdata(a).waveforms))
                            fprintf(awgdata(a).awg, sprintf('WLIS:WAV:NEW "%s",%d,INT', name, npts));
                            awgdata(a).waveforms{end+1}=name;
                            err = query(awgdata(a).awg, 'SYST:ERR?');
                            if ~isempty(strfind(err,'E11113'))
                                fprintf(err(1:end-1));
                                error('Error loading waveform; AWG is out of memory.  Try awgclear(''all''); ');
                            end
                        end
                        awgloadwfm(a,grp.pulses(i).data(dind).wf(j,:), uint16(grp.pulses(i).data(dind).marker(j,:)), name, ch, 0);
                        zerolen{a}(ind(i), j) = -npts;
                        nzpls=1;
                    else
                        zerolen{a}(ind(i), j) = npts;
                    end
                end
            end
            % If no non-zero pulses were loaded, make a dummy waveform so awgclear
            % knows this group was in memory.
            if nzpls == 0
                name=sprintf('%s_1_1',grp.name);
                npts=256;
                if isempty(strmatch(name,awgdata(a).waveforms))
                    fprintf(awgdata(a).awg, sprintf('WLIS:WAV:NEW "%s",%d,INT', name, npts));
                    awgdata(a).waveforms{end+1}=name;
                end
                awgloadwfm(a,zeros(1,npts), zeros(1,npts), name, 1, 0);
            end
            % if the pulse group is added, update it's load time.
            ind=awggrpind(grp.name);
            if ~isnan(ind)
                for i=1:length(awgdata)
                    awgdata(i).pulsegroups(ind).lastload=now;
                end
            end
            
            if dosave
                awgsavedata;
            end   
        end
        
        function load_raw_wf(awg, data, marker, name, chan,define)
            % Send waveform 'data,marker' to the awg with name 'name' intended for channel c.
            % data is scaled and offset by awgdata.scale and awgdata.offset *before* sending.
            start = now;
            if exist('define','var') && define
                fprintf(awg.inst, sprintf('WLIS:WAV:NEW "%s",%d,INT', name, length(data)));
                awg.waveforms{end+1}=name;
            end
            chunksize=65536;
            if(size(data,1) > size(data,2))
                data=data';
            end
            data=(awgdata(a).offset(min(chan,end)) + data)./awgdata(a).scale(chan) + 1;
            tb=find(data > 2);
            tl=find(data < 0);
            if ~isempty(tb) || ~isempty(tl)
                fprintf('Pulse "%s", channel %i, exceeds allowed range: %g - %g (0-2 allowed)\n',name,chan,min(data),max(data));
                data(tb) = 2;
                data(tl) = 0;
            end 
            data = uint16(min(data*(2^(awg.bits-1) - 1), 2^(awg.bits)-1)) + uint16(marker) * 2^awg.bits;
            npts = length(data);
            for os=0:chunksize:npts
                if os + chunksize >= npts
                    fwrite(awg.inst, [sprintf('WLIS:WAV:DATA "%s",%d,%d,#7%07d', name, os, npts-os,2 * (npts-os)),...
                        typecast(data((os+1):end), 'uint8')]);
                else
                    fwrite(awg.inst, [sprintf('WLIS:WAV:DATA "%s",%d,%d,#7%07d', name, os, chunksize,2 * chunksize),...
                        typecast(data((os+1):(os+chunksize)), 'uint8')]);
                end
                fprintf(awg.inst,'');
            end
            time=(now-start)*24*60*60;
            fprintf('Load time: %g seconds for %g points (%g bytes/sec)\n',time, npts,npts*2/time);
        end
        
        function add_trig_pls(awg)
           if awg.slave
               return
           end
           ind = 1+size(awg.waveforms,1);
           for ch = 1:length(awg.chans)
               if ch == awg.trig_pls.chan
                   awg.waveforms{ind,ch} = awg.trig_pls.name;
                   awg.pls_lens(ind,ch) = awg.trig_pls.len;
               else
                   awg.add_zero_pls(ch,awg.trig_pls.len,ind);
               end
           end
        end
        
        function add_zero_pls(awg,chan,len,grp_ind)
           awg.waveforms{grp_ind,chan} = sprintf(awg.zero_pls.zero_name,awg.zero_pls.zero_chan(chan),len);
            return; %FIXME 
        end
        
        function write_seq_file(awg)
            if isempty(awg.waveforms), return, end
            awg.stop();
            empties = cellfun(@isempty,(awg.waveforms));
            for ch = 1:size(empties,2)
               for ll = 1:size(empties,1)
                  if empties(ll,ch)
                      chan_ind = find(~empties(ll,:),1);
                      awg.waveforms{ll,ch} = sprintf('zero_%d_%04d',awg.zero_pls.zero_chan(ch),awg.pls_lens(ll,chan_ind));
                  end
               end
            end
            
            nchan = size(awg.waveforms,2);
            nlines = sum([awg.pulsegroups.n_lines]);
            
            buf =  sprintf('MAGIC 300%1d\r\nLINES %04d\r\n', nchan, nlines); 
            
            for g = 1:length(awg.pulsegroups)
                for i = 1:(awg.pulsegroups(g).n_lines)
                    if i==1 && awg.pulsegroups(g).use_trig && ~awg.slave
                         filename_char = [];
                        for ch=1:nchan
                            filename_char = sprintf('%s"%s",',filename_char,awg.trig_pls.name);
                            reps = 1;
                            jump = awg.pulsegroups(g).st_ind+1;
                            trig = 0;
                        end
                    else
                        filename_char = [];
                        for ch=1:nchan
                            filename_char = sprintf('%s"%s",',filename_char,...
                                awg.waveforms{i+awg.pulsegroups(g).st_ind-1,ch});
                        end
                        % go forward one line only if triggering (i.e. no
                        % trig pls if this is a slave awg)
                        reps=awg.pulsegroups(g).n_rep(i-(awg.pulsegroups(g).use_trig*(~awg.slave)));
                        jump=awg.pulsegroups(g).jump(i-(awg.pulsegroups(g).use_trig*(~awg.slave)));
                        if i ==1 && awg.slave
                            trig = 1;
                        else
                            trig = 0;
                        end
                        %         line = sprintf('"%s","%s",%d,%d,0,1\r\n',...
                        %             ch1_wfm,ch2_wfm,nreps,bTrig);

                    end
                    if isinf(reps)
                        reps = 0; 
                    end
                    line = sprintf('%s%d,%d,%d,%d\r\n',...
                        filename_char,reps,trig,0,0); %0 is logic jump
                        %filename_char,reps,trig,3,jump); %0 is logic jump
                    buf = [buf line]; %#ok<AGROW>
                    
                end
            end
            %buf = [buf, sprintf('TABLE_JUMP 0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1\r\n')]; %dummy event stuff
            %buf = [buf, sprintf('LOGIC_JUMP 1,-1,-1,-1\r\n')];
            %buf = [buf, sprintf('JUMP_MODE SOFTWARE\r\nJUMP_TIMING ASYNC\r\nSTROBE 0\r\n')];
            buf = [buf, sprintf('JUMP_MODE SOFTWARE\r\n')];
            
            fid = fopen([awg.wf_dir,'0000.seq'],'Wb');
            fwrite(fid,buf);
            fclose(fid);
        end
        
        function update_offset(awg,offst,chans)
           if ~exist('chans','var') || isempty(chans)
               chans = 1:length(awg.chans);
           end
           if length(chans)~=length(offst)
              error('size of channel %d different from size of offset %d',size(chans),size(offst)); 
           end
           awg.offst(chans) = offst;
           [ofs,~, of_inds] = unique(awg.offst);
           if any(awg.zero_pls.zero_chan ~=of_inds)
               fprintf('changing offsets. need to remake groups \n');
           end
           awg.zero_pls.zero_chan = of_inds;
        end
        
        function out = load_seq(awg,seq_file)
            if isempty(awg.inst)
               fprintf('empty inst. nothing will be loaded'); 
               out = 0;
               return
            end
           if ~exist('seq_file','var')|| isempty(seq_file)
              seq_file = [awg.wf_dir,'0000.seq']; 
           end
           if isempty(strfind(seq_file,awg.wf_dir))
              seq_file = [awg.wf_dir,seq_file]; 
           end
           seq_file = strrep(seq_file,'/','\\'); %awg likes double slashes
           %give it channel 1. doesn't really matter
           msg = sprintf('SOUR%d:FUNC:USER "%s"',1,seq_file);
           fprintf('%s\n',msg);
           fprintf(awg.inst,msg);
           %fprintf(awg.inst,'SOUR%d:FUNC:USER "z:\\Mikey\\software_dev\\test\\tmp\\awg1\\0000.seq"',1);
           %awg.wait()
           %1ns wait hard coded for slaved awgs
           if awg.slave
              for g = 1:length(awg.pulsegroups)
                 fprintf(awg.inst, sprintf('SEQ:ELEM%d:TWAIT 1\n', awg.pulsegroups(g).st_ind)); 
              end
           end
           awg.write_gotos;
           out = 1;
           
        end
        
        function write_gotos(awg)
            for g = 1:length(awg.pulsegroups)
                inds = awg.pulsegroups(g).st_ind+(0:awg.pulsegroups(g).n_lines-1);
                if awg.pulsegroups(g).use_trig && ~awg.slave
                    gotos = [awg.pulsegroups(g).jump(end), awg.pulsegroups(g).jump];
                elseif awg.slave
                    gotos = awg.pulsegroups(g).jump;
                else
                    gotos = [awg.pulsegroups(g).jump(1), awg.pulsegroups(g).jump];
                end
                for l = 1:length(inds)
                    fprintf(awg.inst,sprintf('SEQ:ELEM%d:GOTO:STAT 0',inds(l)));
                    fprintf(awg.inst,sprintf('SEQ:ELEM%d:GOTO:IND %d',inds(l),gotos(l)));
                    fprintf(awg.inst,sprintf('SEQ:ELEM%d:GOTO:STAT 1',inds(l)));
                end
            end
        end
        
        function on(awg,chans)
            if ~exist('chans','var')|| isempty(chans)
                chans = awg.chans;
            end
           for i = chans
               fprintf(awg.inst, 'OUTPUT%i:STAT 1', i);
           end 
        end
        
        function off(awg,chans)
             if ~exist('chans','var')|| isempty(chans)
                chans = awg.chans;
            end
            for i = chans
                fprintf(awg.inst, 'OUTPUT%i:STAT 0', i);
            end
        end
        
        function stop(awg)
            fprintf(awg.inst, 'AWGC:STOP'); 
        end
        
        function start(awg)
           fprintf(awg.inst, 'AWGC:RUN');
           awg.wait();
        end
        
        function wait(awg) %set the timeout to 5min, query for wait, reset timeout
           to = awg.inst.timeout;
           awg.inst.timeout = 600;
           query(awg.inst, '*OPC?');
           awg.inst.timeout = to; 
        end
        
        function raw(awg,chans)
            if any(awg.israw())
                if isempty(chans)
                    chans = awg.chans;
                end
                    if strcmp(awg.type,'5k')
                        for i = chans
                            fprintf(awgdata(a).awg, 'AWGC:DOUT%i:STAT 1', i);
                        end
                    end
            else
                fprintf('Already raw\n');
            end
        end
        
        function amp(awg,chans)
            if 1% FIXME any(awg.israw())
                if ~exist('chans','var')|| isempty(chans)
                    chans = awg.chans;
                end
                if (awg.is5k)
                    for i = chans
                        fprintf(awg.inst, 'AWGC:DOUT%i:STAT 0', i);
                    end
                end
            else
                fprintf('Already raw\n');
            end
        end
        
        function val = israw(awg,chans)
           val=[];
           if ~exist('chans','var')|| isempty(chans)
               chans = awg.chans;
           end
                if (awg.is5k)
                    for i = chans
                        fprintf(awg.inst, 'AWGC:DOUT%i:STAT?',i);
                        val(end+1) = fscanf(awg.inst,'%f');
                    end
                end  
        end
        
        function val = ison(awg) %if instrument waiting for trigger, returns .5
           val=[];                                
           fprintf(awg.inst, 'AWGC:RST?');
           val(end+1) = .5*fscanf(awg.inst,'%f');                       
        end
        
        function exton(awg,chans)
            if isempty(chans)
                chans = awg.chans;
            end
           for i = chans
               fprintf(awgdata(a).awg, 'SOUR%i:COMB:FEED "ESIG"', i);
           end 
        end
        
        function extoff(awg,chans)
            if isempty(chans)
                chans = awg.chans;
            end
            for i = chans
                fprintf(awgdata(a).awg, 'SOUR%i:COMB:FEED ""', i);
            end
        end
        
        function err(awg)
            er=query(awg.inst, 'SYST:ERR?');
            if strcmp(er(1:end-1), '0,"No error"')
                % Supress blank error messages.
            else
                fprintf(': %s\n',er);
            end
        end
        
        function clr(awg)
            i = 0;
            err2 = sprintf('n/a.\n');
            while 1
                err = query(awg.inst, 'SYST:ERR?');
                if strcmp(err(1:end-1), '0,"No error"')
                    if i > 0
                        fprintf('awg: %i errors. Last %s', i, err2);
                    end
                    break;
                end
                err2 = err;
                i = i + 1;
            end
        end
        
    end
    
end


