%%
s.name = 'testpg';
s.pulses = 2;
s.params = {[4 .2 .3]};
s.varpar{1} =struct('param_num',2,'val',1:10);
%s.dict = {'right'};
s.dict = {struct('prep','@adprep','read','@adread'),'right'};
%%
clear pg
pg = sm_pulsegroup(s);

%% make a new plsdata
global plsdata
plsdata2 = plsdata;
plsdata2.pulses = plsdata2.pulses(14:16);
for j = 1:length(plsdata2.pulses)
    pp=plsdata2.pulses(j).pardef;
    pardef = struct();
    for kk = 1:size(pp,1)
        pardef(kk).elem_num = pp(kk,1);
        if pp(kk,2)<0
            pardef(kk).par = 'time';
        else
            pardef(kk).par = 'val';
        end
        pardef(kk).ind = abs(pp(kk,2));
    end
    plsdata2.pulses(j).pardef = pardef;
    plsdata2.pulses(j).trafofn = func2str(plsdata2.pulses(j).trafofn);
    plsdata2.pulses(j).fill_elem = [];
end
%plsdata2.pulses = rmfield(plsdata2.pulses,'format');
clear plsdata;
global plsdata;
plsdata = plsdata2;
%%
clear s;
s.groups = {pg,pg2};
cg = sm_combo_group(s);

%%
clear s
s.wf_dir= '/Users/michaelshulman/Documents/MATLABL/tmp/';
awg = sm_awg(s);
%%
chans = struct();
[chans(1:4).awg]=deal(1);
[chans(1:4).scale]=deal(10);
[chans(1:4).offset]=deal(0);
for j =1:length(chans)
   chans(j).channel = j; 
end
clear s;
s.awg = awg;
s.chans = chans;
s.name = 'awgdata';
global awgdata;
awgdata = sm_awg_data(s);

%%
awgdata.queue{1} = pg;
awgdata.queue{2} = cg;
%%
clear all;
clear all global
b = instrfind;
delete(b(:));
%plsdata2=load('plsdata_2012_08_22'); 
plsdata2 = load('z:/qDots/awg_pulses/plsdata_2012_08_22.mat');
plsdata2 = plsdata2.plsdata;
if 1
plsdata2.pulses = plsdata2.pulses(14:16);
for j = 1:length(plsdata2.pulses)
    pp=plsdata2.pulses(j).pardef;
    pardef = struct();
    for kk = 1:size(pp,1)
        pardef(kk).elem_num = pp(kk,1);
        if pp(kk,2)<0
            pardef(kk).par = 'time';
        else
            pardef(kk).par = 'val';
        end
        pardef(kk).ind = abs(pp(kk,2));
    end
    plsdata2.pulses(j).pardef = pardef;
    plsdata2.pulses(j).trafofn = func2str(plsdata2.pulses(j).trafofn);
    plsdata2.pulses(j).fill_elem = [];
end
plsdata2.pulses = rmfield(plsdata2.pulses,'format');
clear plsdata; plsdata = plsdata2;

% for j = (1:3)
%     s= struct();
%     s.data = plsdata2.pulses(13+j).data;
%     s.name = plsdata2.pulses(j+13).name;
%     s.pardef = plsdata2.pulses(j+13).pardef;
%     s.trafofn = func2str(plsdata.pulses(13+j).trafofn);
%     plsdata.pulses(j) = sm_pulse(s);
% end
%plsdata = plsdata2;
else
for j = 1:3
   plsdata.pulses(j) = sm_pulse(plsdata2.pulses(13+j));
end
global plsdata
plsdata.tbase = plsdata2.tbase;
end
%plsdata.grpdir = '/Users/michaelshulman/Documents/MATLAB/';
%plsdata.grpdir = '/Users/mikey/Documents/MATLAB/';
plsdata.grpdir = 'z:/Mikey/software_dev/test/';
global plsdata

clear s;
s.name = 'testpg';
s.pulses = 2;
s.params = {[4 .2 .3]};
s.chan = [3 4];
s.varpar{1} =struct('param_num',2,'val',(1:100)*1e-3);
%s.dict = {'right'};
s.dict = {struct('prep','@adprep','read','@adread'),'right'};
clear pg
pg = sm_pulsegroup(s);

s.chan = [4 5];
s.params{1}(2) = .4;
pg2 = sm_pulsegroup(s);

s.params{1}(2) = .23;
s.chan = [5 6];
pg3 = sm_pulsegroup(s);

clear s;
s.groups = {pg,pg2};
cg = sm_combo_group(s);

clear s
s.wf_dir= [plsdata.grpdir,'tmp/awg1/'];
s.trig_pls.name = 'trig_awg1';
s.trig_pls.len = 1000;
s.trig_pls.chan = 2;
s.trig_pls.marker = 1;
awg = sm_awg(s);
s.wf_dir = [plsdata.grpdir, 'tmp/awg2/'];
s.trig_pls.name = 'trig_awg2';
s.slave = 1;
awg2 = sm_awg(s);

chans = struct();
for j =1:8
    if j <5
        chans(j).awg = 1;
       chans(j).channel = j; 
    else
        chans(j).awg = 2;
        chans(j).channel = j-4;
    end
    chans(j).scale =10;
    chans(j).offset =0;
end

clear s;
s.awg = awg;
s.awg(2) = awg2;
s.chans = chans;
s.name = 'awgdata';
global awgdata;
awgdata = sm_awg_data(s);


awgdata.queue{1} = pg;
awgdata.queue{2} = pg2;%cg;
awgdata.queue{3} = pg3;
try
    awgdata.awg(1).inst= visa('ni','tcpip::140.247.189.8::INSTR');
    fopen(awgdata.awg(1).inst);
end

%%
clear all;
clear all global
b = instrfind;
delete(b(:));
%plsdata2=load('plsdata_2012_08_22'); 
plsdata2 = load('z:/qDots/awg_pulses/plsdata_2012_08_22.mat');
plsdata2 = plsdata2.plsdata;
if 1
plsdata2.pulses = plsdata2.pulses(14:end);
for j = 1:length(plsdata2.pulses)
    pp=plsdata2.pulses(j).pardef;
    pardef = struct();
    for kk = 1:size(pp,1)
        pardef(kk).elem_num = pp(kk,1);
        if pp(kk,2)<0
            pardef(kk).par = 'time';
        else
            pardef(kk).par = 'val';
        end
        pardef(kk).ind = abs(pp(kk,2));
    end
    plsdata2.pulses(j).pardef = pardef;
    if ~isempty(plsdata2.pulses(j).trafofn)
        plsdata2.pulses(j).trafofn = func2str(plsdata2.pulses(j).trafofn);
    else
        plsdata2.pulses(j).trafofn = (plsdata2.pulses(j).trafofn);
    end
    plsdata2.pulses(j).fill_elem = [];
end
plsdata2.pulses = rmfield(plsdata2.pulses,'format');
clear plsdata; plsdata = plsdata2;

% for j = (1:3)
%     s= struct();
%     s.data = plsdata2.pulses(13+j).data;
%     s.name = plsdata2.pulses(j+13).name;
%     s.pardef = plsdata2.pulses(j+13).pardef;
%     s.trafofn = func2str(plsdata.pulses(13+j).trafofn);
%     plsdata.pulses(j) = sm_pulse(s);
% end
%plsdata = plsdata2;
else
for j = 1:3
   plsdata.pulses(j) = sm_pulse(plsdata2.pulses(13+j));
end
global plsdata
plsdata.tbase = plsdata2.tbase;
end
%plsdata.grpdir = '/Users/michaelshulman/Documents/MATLAB/';
%plsdata.grpdir = '/Users/mikey/Documents/MATLAB/';
plsdata.grpdir = 'z:/Mikey/software_dev/test/';
global plsdata

clear s;
s.name = 'testpg';
s.pulses = 2;
s.params = {[4 .2 .3]};
s.chan = [3 4];
s.varpar{1} =struct('param_num',2,'val',(1:100)*1e-3);
%s.dict = {'right'};
s.dict = {struct('prep','@adprep','read','@adread'),'right'};
clear pg
pg = sm_pulsegroup(s);

% s.chan = [4 5];
% s.params{1}(2) = .4;
% pg2 = sm_pulsegroup(s);
% 
% s.params{1}(2) = .23;
% s.chan = [5 6];
% pg3 = sm_pulsegroup(s);

% clear s;
% s.groups = {pg,pg2};
% cg = sm_combo_group(s);



clear s
s.wf_dir= [plsdata.grpdir,'tmp/awg1/'];
s.trig_pls.name = 'trig_awg1';
s.trig_pls.len = 1000;
s.trig_pls.chan = 2;
s.trig_pls.marker = 1;
awg = sm_awg(s);
%s.wf_dir = [plsdata.grpdir, 'tmp/awg2/'];
%s.trig_pls.name = 'trig_awg2';
%s.slave = 1;
%awg2 = sm_awg(s);

chans = struct();
for j =1:4
    if j <5
        chans(j).awg = 1;
       chans(j).channel = j; 
    else
        chans(j).awg = 2;
        chans(j).channel = j-4;
    end
    chans(j).scale =10;
    chans(j).offset =0;
end

clear s;
s.awg = awg;
%s.awg(2) = awg2;
s.chans = chans;
s.name = 'awgdata';
global awgdata;
awgdata = sm_awg_data(s);


awgdata.queue{1} = pg;


plschans = {[2 1], [3 4]};
for j = 1:length(plschans)
%plen, max_evo, eps, evo
for eps = linspace(1,4,10)
clear s;
s.options = 'pack loop';
s.name = 'ramsey';
s.pulses = [57 57];
s.params = {[4, .11, eps 0],[4, .11, eps+.1 0]};
s.chan = plschans{j};
s.varpar{1} =struct('param_num',4,'val',(1:100)*1e-3);
s.varpar{2} =struct('param_num',4,'val',(1:100)*1e-3);
%s.dict = {'right'};
s.dict = {struct('prep','@adprep','read','@adread'),'right'};
clear pg
awgdata.queue{end+1} = sm_pulsegroup(s);
end
end

try
    awgdata.awg(1).inst= visa('ni','tcpip::140.247.189.8::INSTR');
    fopen(awgdata.awg(1).inst);
end