%%
clear all
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
s.elems = [st,rl,ap,sep,ar,wt,rd];
s.pardef(1).elem_num = 4; s.pardef(1).par = 'time'; s.pardef(1).ind = 1;
s.pardef(2) = s.pardef(1); s.pardef(2).par = 'mult';
s.fill.elem = 6; s.fill.time = 4;
s.trafofn = '@(x)[x(1)*5e-3, x(2)]';
pls = pls_pls(s);


%%

pg.name = 'test';
pg.pulses = pls;
pg.params = [.1,1.2];
pg.varpar.param_num = 1;
pg.varpar.val = 1:100;
pg.chan = [3,4];
pgroup = pls_plsgroup(pg);


