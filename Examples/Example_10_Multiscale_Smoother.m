% Example_10_Multiscale_Smoother: 
% An idealized test for multiscale nesting using boxes with a large min_el
% ratio

bbox = [0, 1; 0,1]; 
boubox = [0,0; 1,0; 1,1; 0,1; 0,0; NaN NaN ];
min_el = 1e3; 

gdat1 = geodata('pslg',boubox,'bbox',bbox,'h0',min_el);
fh1 = edgefx('geodata',gdat1,'g',0.2);

bbox2 = [-1, 2; -1,2]; 
boubox2 = [-1,-1; 2,-1; 2,2; -1,2; -1,-1; NaN NaN ];
min_el2 = min_el*10; 

gdat2 = geodata('pslg',boubox2,'bbox',bbox2,'h0',min_el2);
fh2 = edgefx('geodata',gdat2,'g',0.2);

mshopts = meshgen('ef',{fh2, fh1},'bou',{gdat2,gdat1},...
                  'plot_on',1,'qual_tol',0.0025,'cleanup',0);
mshopts = mshopts.build;

m = mshopts.grd; 

plot(m)

m = m.clean;

plot(m)