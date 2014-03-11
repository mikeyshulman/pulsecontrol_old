classdef sm_awg_data < handle
    %sm_awg_data < handle, wrapper class for admining multiple awgs
    %   has properties: 
        %name: mostly useless
        % awg: array of type sm_awg
        % chans: struct array of human readable channles with fields:
            % awg (which awg in awgdata.awg) 
            % channel: which chennel in that awg
            % scale: full scale of channel, in mV
            %: offset: in mV
        % memory: struct array of groups loaded in memory with fields:
            % grp: of type sm_plsgrp or the like
            % awg_ind: which awg(s) is it on
            % st_ind: same size as awg_ind, where seq starts on each awg
            % n_lines: number of seq lines it takes on each awg
        % queue: a cell array of type sm_pulsegrp or the like
            % these will get dumped to memory when awgdata.write_seq is
            % called
          
        % Some important methods:
            % constructor takes a struct with the same fields
            % cntrl('str'): like old awgcntrl, for start/stopping awgs
                %e.g. awgdata.cntrl('on start wait err')
            % write_seq: write the queue to a seq file
            % lint: check things for consistency
            % load_seq: load the sequence onto each awg
            
    properties
        name='awgdata';
        awg;
        chans = struct('awg',[],'channel',[],'scale',[],'offset',[]);
        memory = struct('grp',{},'awg_ind',{},'seq_ind',{},'n_lines',{});
        queue = {};
    end
    
    properties (Transient = true)
       changed = 0; 
    end
    
    methods
        function ad = sm_awg_data(s)
           fn = fieldnames(s);
           for j = 1:length(fn)
              ad.(fn{j}) = s.(fn{j}); 
           end
        end
        
        function out = cntrl(awgdata,ctl_str)
            out = [];
            c2=regexp(ctl_str,'[^\t][^ \t]*','match');
            for j = 1:length(c2)
               for k = 1:length(awgdata.awg)
                  awgdata.awg(k).(c2{j});
                   %out = [out, awgdata.awg.(c2{j})()];  %#ok<AGROW>
               end
            end
        end
        
        function sync(awgdata)
            awgdata.write_seq();
            awgdata.load_seq();
        end
        
        function status = lint(awgdata)
           [~,unique_inds,ii] = unique({awgdata.queue.name});
           if length(unique_inds) < length(ii)
              error('repeated names of groups %i',instersect(unique_inds,ii)); 
           else
               status = 1;
           end
           
        end
        
        function update_scale(awgdata,new_scale,ch)
           if ~exist('ch','var')|| isempty(ch)
               ch = 1:length(awgdata.chans);
           end
           if length(ch)~=length(new_scale)
              error('length of channel not same as new scale'); 
           end
           old_scale = [awgdata.chans(ch).scale];
           if any(old_scale~=new_scale)
               fprintf('scales changed!!! remake wfs\n');
               awgdata.changed = 1;
           end
           for c = 1:length(ch)
              awgdata.chans(ch(c)).scale = new_scale(c);
           end
        end
        
        function update_offset(awgdata,new_off,ch)
            if ~exist('ch','var')|| isempty(ch)
                ch = 1:length(awgdata.chans);
            end
            if length(ch)~=length(new_off)
                error('length of channel not same as new offset');
            end
            old_off = [awgdata.chans(ch).offset];
            if any(old_off~=new_off)
                awgdata.changed = 1;
                fprintf('offsets changed. remake wfs! \n');
            end
            for c = 1:length(ch)
               awgdata.awg(awgdata.chans(c).awg).update_offset(new_off(c),ch(c));
               awgdata.chans(c).offset = new_off(c);
            end            
        end
        
        function write_seq(awgdata,opts)
            global plsdata;
            if ~exist('opts','var')
                opts = '';
            end
            %awgdata.lint()
            grp_changed = true(1,length(awgdata.queue));
            if isempty(strfind(opts,'force'))
                for j = 1:length(awgdata.memory)
                    grp_changed(j) = ~isequal(awgdata.queue{j}.to_struct,awgdata.memory(j).grp);
                end
            end
            awgdata.changed = any(grp_changed);
            start_ind = find(grp_changed,1);
            if isempty(start_ind)
                fprintf(['queue is same as memory. nothing to do\n',...
                    'consider ''force'' option\n']);
                return
            end
            fprintf('starting sync from group %i\n',start_ind);
            
            use_trig = zeros(1,length(awgdata.queue));
            for j = 1:length(awgdata.queue)
               use_trig(j) = (isempty(awgdata.queue{j}.n_rep) || ...
                   awgdata.queue{j}.n_rep(1) ~= Inf) && isempty(strfind(awgdata.queue{j}.options, 'notrig')); 
            end
            for a = 1:length(awgdata.awg) %clear stuff out
                if isempty(awgdata.awg(a).pulsegroups)
                    chg = 1;
                else
                    chg = find([awgdata.awg(a).pulsegroups.grp_ind]==start_ind);
                end
                if isempty(chg)
                    chg = 1;
                end
                %awgdata.awg(a).pulsegroups(chg:end)=[];
                if chg ==1
                    awgdata.awg(a).waveforms = cell(0,length(awgdata.awg(a).chans));
                    awgdata.awg(a).pls_lens = zeros(size(awgdata.awg(a).waveforms));
                else
                    awgdata.awg(a).waveforms(awgdata.awg(a).pulsegroups(chg).st_ind:end,:)=[];
                    awgdata.awg(a).pls_lens(awgdata.awg(a).pulsegroups(chg).st_ind:end,:)=[];
                end
            end
            %now add a trigger pls to each awg if there are no wfs
            for a = 1:length(awgdata.awg) %add triggers to each
                if 1 %awgdata.awg(a).is_wf_list_empty
                clear wf_data;
                wf_data.wf=zeros(1,awgdata.awg(a).trig_pls.len);
                wf_data.marker=repmat(1,1,awgdata.awg(a).trig_pls.len); % channel 1&3 marker 1
                %write_wf_file(wf_data,filename,chan_ind,clock_f)
                fname = [awgdata.awg(a).wf_dir,sprintf(awgdata.awg(a).trig_pls.name)];
                write_wf_file(wf_data,fname,1,awgdata.awg(a).clk);
%                awgloadwfm(a,zdata,zmarker,sprintf('trig_%08d',awgdata(a).triglen),1,1);
                %awgdata.awg(j).add_trig_wf();
                end
            end
            go_tos = zeros(0,length(awgdata.awg));
            n_reps = zeros(0,length(awgdata.awg));
            %jumps = zeros(0,length(awgdata.awg));
            basename = 'awg%d_grp%04d_ch%04d_wf%04d.wf';
            seq_lines = zeros(1,length(awgdata.awg));
            awgdata.memory(start_ind:end) = [];
            for g = start_ind:length(awgdata.queue)
                use_awg =false(1,length(awgdata.awg));
                use_awg([awgdata.chans(awgdata.queue{g}.chan).awg])=true;
                if 1%awgdata.changed || isempty(awgdata.queue{g}.pack_data) || ...
                        %isempty(awgdata.queue{g}.pack_data(1).wf)||isempty(awgdata.queue{g}.pack_data(1).marker)
                   %pg=make(plsgrp,ind, clk, tbase, time)
                   clks = [awgdata.awg(use_awg).clk];
                   if any(clks ~=clks(1))
                       error('different awgs with different clocks for same group!');
                   end
                    pg_struct=awgdata.queue{g}.make([],clks(1),plsdata.tbase,[]); 
                end
                if isempty(awgdata.queue{g}.n_rep)
                   nr = ones(1,size(pg_struct.pack_data,2)); 
                else
                    nr = awgdata.queue{g}.n_rep;
                end
                if isempty(awgdata.queue{g}.jump)
                    jmp = [2:size(pg_struct.pack_data,2), 1];
                else
                    jmp = awgdata.queue{g}.jump;
                end
                if use_trig(g)
                    seq_lines = seq_lines + use_awg.*abs([awgdata.awg.slave]-1); %add seq line for trigger
                    if isempty(go_tos)
                       go_tos(1,use_awg) = 1; 
                    else
                        go_tos(end+1,use_awg) = go_tos(end,use_awg)+1; %#ok<AGROW> % got to next line
                    end
                    n_reps(end+1,use_awg) = 1;%#ok<AGROW> 
                end
                awg_ind = 0;
                for ch = (awgdata.queue{g}.chan)
                    if awgdata.chans(ch).awg ~= awg_ind
                        first_time_loaded = 1;
                    else
                        first_time_loaded =0;
                    end
                    awg_ind = awgdata.chans(ch).awg;
                    pack = ~isempty(strfind(awgdata.queue{g}.options,'pack'));
                    loop = ~isempty(strfind(awgdata.queue{g}.options,'loop'));
                    once = ~isempty(strfind(awgdata.queue{g}.options,'once'));
                    
                    if pack
                        if any(diff(jmp(1:end-1))~=1)
                            error('cannot jump within grp if packing wfs');
                        end
                        if any(nr~=1)
                           error('cannot repeat wfs with pack option'); 
                        end
                        pack_data = struct();
                        pack_data.wf = [pg_struct.pack_data.wf];
                        pack_data.marker = [pg_struct.pack_data.marker];
                    else
                        pack_data = pg_struct.pack_data;
                    end
                    if first_time_loaded
                       st_ind = 1+size(awgdata.awg(awg_ind).waveforms,1);
                    end
                    seq_line_tmp = seq_lines;%seq_lines(:,end);
                    if use_trig(g) && first_time_loaded
                       awgdata.awg(awg_ind).add_trig_wf;
                       awgdata.awg(awg_ind).add_trig_pls;
                    end
                    for wfind = 1:length(pack_data)
                        wf_name = sprintf(basename,awg_ind,g,ch,wfind);
                        chan_ind = find(awgdata.queue{g}.chan==ch);
                        pack_data(wfind).wf(chan_ind,:) = (pack_data(wfind).wf(chan_ind,:)+awgdata.chans(ch).offset)./awgdata.chans(ch).scale; %#ok<AGROW>
                        too_high = pack_data(wfind).wf(chan_ind,:)>1;
                        too_low = pack_data(wfind).wf(chan_ind,:) < -1;
                        if any(too_high)|| any(too_low)
                            fprintf('warning, pulse %i of group %i on channel %i out of range. clipping...\n',wfind,g,ch);
                        end
                        pack_data(wfind).wf(chan_ind,too_high) = 1;%#ok<AGROW>
                        pack_data(wfind).wf(chan_ind,too_low) = -1; %#ok<AGROW>% rectify signal
                        write_wf_file(pack_data(wfind),[awgdata.awg(awg_ind).wf_dir,wf_name],chan_ind,awgdata.awg(awg_ind).clk);
                        awgdata.awg(awg_ind).waveforms{-awgdata.awg(awg_ind).slave+st_ind+wfind,awgdata.chans(ch).channel} = wf_name;
                        awgdata.awg(awg_ind).pls_lens(-awgdata.awg(awg_ind).slave+st_ind+wfind,awgdata.chans(ch).channel) = length(pack_data(wfind).wf);
                        if first_time_loaded
                            seq_lines(awg_ind) = seq_lines(awg_ind)+use_awg(awg_ind);% add one seq line for each awg used
                        end
                    end
                    if first_time_loaded
                        for a = awg_ind%unique(find(use_awg))
                            awgdata.awg(a).pulsegroups(end+1).grp = awgdata.queue{g};
                            awgdata.awg(a).pulsegroups(end).st_ind = st_ind;
                            if pack
                                if loop
                                    awgdata.awg(a).pulsegroups(end).jump = st_ind+use_trig(g)*(~awgdata.awg(a).slave);
                                    awgdata.awg(a).pulsegroups(end).n_rep = Inf;
                                else
                                    awgdata.awg(a).pulsegroups(end).jump = st_ind+2;
                                    awgdata.awg(a).pulsegroups(end).n_rep = 1;
                                end
                            else
                                if once
                                   jmp(end) = jmp(end-1)+1; 
                                end
                                awgdata.awg(a).pulsegroups(end).jump = seq_line_tmp(a)+jmp;%seq_line_tmp(min(a,end))+jmp;
                                awgdata.awg(a).pulsegroups(end).n_rep = nr;
                            end
                            awgdata.awg(a).pulsegroups(end).n_lines = use_trig(g)+wfind-awgdata.awg(a).slave;
                            awgdata.awg(a).pulsegroups(end).use_trig = use_trig(g);
                            awgdata.awg(a).pulsegroups(end).grp_ind = g;
                        end
                    end
                end %end c loop
                
                awg_ind = 0;%#ok<NASGU> %reset it so next group will get first_time_loaded
                awgdata.memory(g).grp = awgdata.queue{g}.to_struct();
                awgdata.memory(g).awg_ind = find(use_awg);
                awgdata.memory(g).seq_ind = zeros(size(awgdata.memory(g).awg_ind));
                awgdata.memory(g).n_lines = zeros(size(awgdata.memory(g).awg_ind));
                for a = 1:length(awgdata.memory(g).seq_ind)
                   awgdata.memory(g).seq_ind(a) = awgdata.awg(a).pulsegroups(end).st_ind; 
                   awgdata.memory(g).n_lines(a) = (awgdata.awg(a).pulsegroups(end).n_lines);
                end
            end %end g loop
            write_start = tic;
            for a = 1:length(awgdata.awg)
                %write_seq_file(awg,jumps,reps,trig)
               awgdata.awg(a).write_seq_file();
               fprintf('wrote seq for awg %i / %i. Ellapsed time: %.2f seconds \n',a,length(awgdata.awg),toc(write_start));
            end
            %fprintf('\n');
            if ~isempty(strfind(opts,'purge'))
                awgdata.queue = {};
            end
            
            %now loop through awgs and seq lines and write zero wfs
            % the following is a little silly since there will be some
            % extra zero writing, but it takes negligible time and makes
            % the code nice and readable
            for a = 1:length(awgdata.awg)
                [~,ia,~]=unique(awgdata.awg(a).zero_pls.zero_chan);
               for ll = 1:size(awgdata.awg(a).waveforms,1)
                  for ch = 1:length(awgdata.awg(a).chans)
                     if awgdata.awg(a).pls_lens(ll,ch)==0 && any((ia-ch)==0)% initialized to zero
                         ind = find(awgdata.awg(a).pls_lens(ll,:)>0,1);
                         if ~isempty(ind)
                             wf_data.wf = awgdata.awg(a).offst(awgdata.chans(ch).channel)*ones(1,awgdata.awg(a).pls_lens(ll,ind));
                             wf_data.marker = zeros(size(wf_data.wf));
                             zname = sprintf(awgdata.awg(a).zero_pls.zero_name,ia(ch),length(wf_data.wf));
                             write_wf_file(wf_data,[awgdata.awg(a).wf_dir,zname],1,awgdata.awg(a).clk);
                         end
                     end
                  end
               end
            end
            awgdata.changed = 0;
        end % end function
        
        function success = load_seq(awgdata,seqfile)
            t_start = tic;
            if ~exist('seqfile','var')
                seqfile = cell(1,length(awgdata.awg));
            end
            success = ones(1,length(awgdata.awg));
           for a= 1:length(awgdata.awg)
              success(a) = awgdata.awg(a).load_seq(seqfile{a});
              fprintf('loaded %i of %i seq files in %.2d seconds \n',a,length(awgdata.awg),toc(t_start));
           end
        end
    end
    
end
function out = write_wf_file(wf_data,filename,chan_ind,clock_f)
%write raw waveforms and markers
% wf should be scaled to be in the range [-1 1]
    
    if length(wf_data.wf)<chan_ind 
        out = 0;
        return;
    end
        
%    clock_f = 1e9/wf_data(channel).ClockT; %Hz;
    V = wf_data.wf(chan_ind,:);
    mkrs = wf_data.marker(chan_ind,:);
    
    ndel = mod(length(V),4);
    if ndel
        warning('clipping end of pulses to make npoints divisible by 4');
        V(end-ndel+1:end) = [];
        for i=1:length(mkrs)
            mkrs(i).State(end-ndel+1:end) = [];
        end
    end
    npoints = length(V);
    marker = zeros(1, npoints, 'uint8');
    z_marker = marker;
%     for i=1:length(mkrs)
%         cur_mkr = cast(mkrs(i).State,'uint8');
%         marker = marker + cur_mkr*2^(i-1);
%     end
    cur_mkr = cast(mkrs,'uint8'); marker = mkrs;
%    marker = marker + cur_mkr*2;
    buf = reshape([reshape(typecast(single(V), 'uint8'), 4, npoints); marker],...
            5*npoints, 1);
    
      
    fid = fopen(filename,'Wb', 'ieee-le');
    
    fwrite(fid, [sprintf('MAGIC 1000\r\n#7%07d', 5*npoints), buf', ...
                sprintf('CLOCK %.8e \r\n', clock_f)]);
    fclose(fid);
if 0
    z_marker = zeros(1,npoints, 'uint8');
    z_V = zeros(size(V));
    
    z_buf = reshape([reshape(typecast(single(z_V), 'uint8'), 4, npoints); z_marker],...
            5*npoints, 1);
    [token,remain] = strtok(filename,'.');
    z_filename = [token, 'z', remain];
        
    fid = fopen(z_filename,'Wb', 'ieee-le');
    fwrite(fid, [sprintf('MAGIC 1000\r\n#7%07d', 5*npoints), z_buf', ...
                sprintf('CLOCK %.8e \r\n', clock_f)]);
    fclose(fid); 
end
    out = 1;
end


