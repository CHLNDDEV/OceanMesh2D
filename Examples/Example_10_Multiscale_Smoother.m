% Smooth outer testing...

bbox = [0, 1; 0,1]; 
boubox = [0,0; 1,0; 1,1; 0,1; 0,0; NaN NaN ];
min_el = 1e3; 

gdat1 = geodata('pslg',boubox,'bbox',bbox,'h0',min_el);
fh1 = edgefx('geodata',gdat1,'g',0.25);

bbox2 = [-1, 2; -1,2]; 
boubox2 = [-1,-1; 2,-1; 2,2; -1,2; -1,-1; NaN NaN ];
min_el2 = min_el*3; 

gdat2 = geodata('pslg',boubox2,'bbox',bbox2,'h0',min_el2);
fh2 = edgefx('geodata',gdat2,'g',0.25);

mshopts = meshgen('ef',{fh2, fh1},'bou',{gdat2,gdat1},'plot_on',1,'cleanup',0);
mshopts = mshopts.build;

m = mshopts.grd; 

plot(m)

m = m.clean;

plot(m)