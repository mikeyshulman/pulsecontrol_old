%%
%clear all
rl = pls_reload;
rd = pls_readout;
ap = pls_adprep;
ar = pls_adread;
wt = pls_wait;
st = pls_raw;

rl.data.time=.4;
rl.data.pos = [-1 -2];
rl.data.after_wait = .1;
rl.data.ramp_time = .02;

rd.data.time = 1;
rd.data.ed_dly=0;
rd.data.meas_pt=[0 0];
rd.data.st_dly = .03;

ap.data.time = .3;
ap.data.start = 1;
ap.data.end = 4;
ap.data.mult = [1 -1];

ar.data = ap.data;
ar.data.start = ap.data.end;
ar.data.end = ap.data.start;

wt.data.time = .5;
wt.data.val = [0,0];

sep = pls_wait;
sep.data.val = [4 -4];

st.data.time = 0;
st.data.val = [0;0];

%%

pe = rl;
pe(3) = ap;
pe(5) = rd;

%%

s= struct();
%s.elems = [st,rl,ap,sep,ar,wt,rd];
s.elems = {pls.raw,pls.reload,pls.adprep,pls.wait,pls.adread,pls.wait,pls.readout};
s.pardef(1).elem_num = 4; s.pardef(1).par = 'time'; s.pardef(1).ind = 1;
s.pardef(2) = s.pardef(1); s.pardef(2).par = 'mult';
s.fill.elem = 6; s.fill.time = 4;
s.trafofn = '@(x)[x(1)*5e-3, x(2)]';
pulse = pls_pls(s);
pulse.elems(4).data.val = [4,-4];
pulse.elems(1).data.val = [0;0]; pulse.elems(1).data.time = 0;


%%

pg.name = 'test';
pg.pulses = pulse;
pg.params = [.1,1.2];
pg.varpar.param_num = 1;
pg.varpar.val = 1:100;
pg.chan = [3,4];
pgroup = pls_plsgroup(pg);


%%

dd.dbzprep = pls_wait; dd.dbzprep.data.val = [-4,4];dd.dbzprep.data.time = .005;
dd.dbzread = pls_blank('@dbzprep');
dd.dbzpi = copy(dd.dbzprep); dd.dbzpi.data.time = .01;

ap = pls_adprep;
ap.data.time = .3;
ap.data.start = 1;
ap.data.end = 4;
ap.data.mult = [1 -1];
dd.adprep = copy(ap);
dd.adread = pls_adread;
dd.adread.data =dd.adprep.data;
dd.reload = pls_reload();
dd.reload.data = struct('time',.2,'after_wait',.05,'ramp_time',.005,'pos',[-1,-2]);
dd.readout = pls_readout;
dd.readout.data=struct('time',1,'st_dly',.2,'ed_dly',0,'meas_pt',[],'flag',1);
dd.start = st;

dict = pls_dict('right',dd);

