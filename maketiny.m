clearvars; close all; clc; 

boubox = [-73.9647   40.5696
    -73.9650   40.5569
    -73.9334   40.5569
    -73.9342   40.5731
    -73.9644   40.5731
    -73.9647   40.5696] ;

weir{1} = [-73.9537   40.5568    0.0002
    -73.9541   40.5632         0
    -73.9545   40.5704         0
    -73.9546   40.5731         0];

gdat = geodata('bbox',boubox,'h0',25,'outer',boubox,'mainland',boubox,'weirs',weir) ;

plot(gdat) ; 

fh = edgefx('geodata',gdat,'dis',0.15) ;

plot(fh) ; 

mshopts = meshgen('bou',gdat,'ef',fh,'plot_on',1,'proj','equi') ;
mshopts = mshopts.build ; 
m = mshopts.grd ; 

plot(m,'tri') ; 

%% 
m = makens(m,'outer',0) ; 

m = makens(m,'weirs',gdat) ; 
%% 
m.b = m.p(:,1)*0 + 3 ; 

m = Calc_tau0(m) ; 

pos = [-73.9534   40.5566
  -73.9536   40.5589
  -73.9543   40.5716
  -73.9543   40.5734
  -73.9334   40.5733
  -73.9329   40.5563] ; 

% make the right side all dry at the start 
in = inpoly(m.p,pos) ; 

ss(:,1) = m.p(:,1)*0 ; 

ss(in,1) = 1; 

m = Calc_f13_inpoly(m,'surface_submergence_state',[],[],[],ss) ;

% jack up the friction on the right too 
cf(:,1) = m.p(:,1)*0 + 0.025 ; 

cf(in,1) = 0.06; 

m = Calc_f13_inpoly(m,'mannings_n_at_sea_floor',[],[],[],cf) ;

%%
m = renum(m) ; 

write(m,'tiny2') ;