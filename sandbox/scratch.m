%%
%clear all
rl = pls_reload;
rd = pls_readout;
ap = pls_adprep;
ar = pls_adread;
wt = pls_wait;
st = pls_raw;


rl.time=.4;
rl.pos = [-1 -2];
rl.after_wait = .1;
rl.ramp_time = .02;

rd.time = 1;
rd.ed_dly=0;
rd.meas_pt=[0 0];
rd.st_dly = .03;

ap.time = .3;
ap.start = 1;
ap.finish = 4;
ap.mult = [1 -1];

ar.time = .3;
ar.start = 4;
ar.finish = 1;
ar.mult = [1 -1];


wt.time = .5;
wt.val = [0,0];

sep = pls_wait;
sep.val = [4 -4];

st.time = 0;
st.val = [0;0];

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
clear dd;
dd.dbzprep = pls_wait; dd.dbzprep.val = [-4,4];dd.dbzprep.time = .005;
dd.dbzread = pls_blank('@dbzprep');
dd.dbzpi = copy(dd.dbzprep); dd.dbzpi.time = .01;

ap = pls_adprep;
ap.time = .3;
ap.start = 1;
ap.finish = 4;
ap.mult = [1 -1];
dd.adprep = copy(ap);
dd.adread = pls_adread;
dd.adread.start =1; dd.adread.finish = 4; dd.adread.time = .4;
dd.reload = pls_reload();
dd.reload.time = .2; dd.reload.after_wait = .05; dd.reload.ramp_time = .005; dd.reload.pos = [-1,-2];
dd.readout = pls_readout;
dd.readout.time = 1; dd.readout.st_dly = .2; dd.readout.ed_dly = 0; dd.readout.meas_pt = [0 0]; dd.readout.flag = 1;
dd.start = st;
dd.exch = pls_wait;


dict = pls_dict('left',dd);

