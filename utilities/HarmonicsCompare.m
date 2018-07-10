function HarmonicsCompare

%-------------------------------------------------------------------------%
%                        Tidal Harmonics Comparison                       %
%                                                                         %
%   This function is designed to create figures for model-station         %
%   harmonics comparison.  It takes as input a text file with the         %
%   following format.                                                     %
%                                                                         %
%   #-of-stations #-of-constituents                                       %
%   'A list of Markers to be used, should equal #-of-constituents'        %
%   station-id  longitude   latitude                                      %
%   const-name  obs-amp     mod-amp     obs-phs     mod-phs               %
%     .                                                                   %
%     .                                                                   %
%     .                                                                   %
%   station-id  longitude   latitude                                      %
%   const-name  obs-amp     mod-amp     obs-phs     mod-phs               %
%     .                                                                   %
%     .                                                                   %
%     .                                                                   %
%                                                                         %
%   OPTIONS:                                                              %
%   runtitle = ___  This string value will be added to the plot titles    %
%                   of the figures.                                       %
%   filetitle = __  This string value will be used to label each of the   %
%                   file produced by the code.                            %
%   flag = _______  This value identifies the flag for any missing        %
%                   data points                                           %
%   sflag = ______  If this value is 1 then the code will produce a text  %
%                   file with a table of the statistical data             %
%   stsused = ____  This is a 6 point vector flagging specific statistics %
%                   i1 - R-Squared                                        %
%                   i2 - Slope of best fit line y=mx                      %
%                   i3 - Standard Deviation of difference in obs-mod      %
%                   i4 - Average difference between obs and mod           %
%                   i5 - Average of the absolute difference               %
%                   i6 - Normalized RMSE                                  %
%   indiv = ______  If this value is 1 then the code will produce figures %
%                   for each individual constituent.  If this value is 0  %
%                   the code will skip individual constituent figures.    %
%   plotall = ____  If this value is 1 then the code will produce a       %
%                   figure will all constituents together.                %
%   bins = _______  This option allows for the separation of the data into%
%                   specific bins, such as regions or states.             %
%   filetypes = __  This is cell input of the filetypes for the figures,  %
%                   they must be the proper string notation for the       %
%                   'print' command in Matlab.                            %
%   xml = ________  If this value is 1 then the code will produce a text  %
%                   file that can be converted into an xml file for use   %
%                   with Google Maps.                                     %
%   ms = _________  Marker size                                           %
%   ampaxis = ____  This is a vector with as many rows as constituents,   %
%                   the first column specifies the desired minimum axis   %
%                   value, the second column specifies the desired max    %
%                   axis value.  This vector is for the amplitudes.       %
%   phsaxis = ____  Same as ampaxis, but for phases.                      %
%   fid = ________  Name of file containing Harmonics results.            %
%   fsz = ________  Fontsize for figures                                  %
%-------------------------------------------------------------------------%
%   Updated by Aaron Donahue to accept String Station ID's  03-11-2012    %
%   Updated by Aaron Donahue to accept bins from a file and               %
%                            to use the 'print' command     04-27-2012    %
%   Updated by Aaron Donahue to fix "All Plots" markers bug 05-04-2012    %
%   Updated by Aaron Donahue to make nicer figures          10-17-2012    %
%   Updated by Aaron Donahue to fix statistics              02-27-2013    %
%   Updated by Keith 2017-
%-------------------------------------------------------------------------%
clear StnId AmpData PhsData
global flag sflag ststable stsused bins mycolors ms myMarks filetypes
global stsvari cntrl AmpData temp2 PhsData StnId errfid
global ampobs ampmod phsobs phsmod mybins Const temp1
global allamp allphs markervec indiv fsz fszstat
global T_o 

%% User Defined Values for Flags
runtitle = 'September 1-November 31 2012';
filetitle = '';
flag = -99999; %note this flag value 
sflag = 1;
stsused = [1 1 1 1 1 1];
indiv   = 0;
plotall =  1;
plotbins = 0;
bins = [];
filetypes = {'-dpng'};
xml = 1 ;
ms = 6;            % Marker Size
fsz = 8;
fszstat = 8;        % Fontsize for stats on figs
ampaxis = [0 0.2;0 0.2;0 0.08;0 0.55;0 0.22;0 0.08;0 0.05;0 0.08];
phsaxis = [-20 200;-120 360;-120 360;-80 360;-120 360;-80 360;-40 360;0 360];
ampext = [min(ampaxis(:,1)) max(ampaxis(:,2))];
phsext = [min(phsaxis(:,1)) max(phsaxis(:,2))];
%% Read the dataset
% fid = fopen('Harm_SL18tx.txt');              % Dataset file
fid = fopen('Harm.txt') ;
T_o = readtable('TIDES_Base_V4_Stations.csv'); 
temp = textscan(fid,'%d',2);
cntrl = temp{1};
temp = textscan(fid,'%s',cntrl(2));
myMarks = temp{1};
% StnId = zeros(cntrl(1),1);
LatLon = zeros(cntrl(1),2);
AmpData = cell(cntrl(1),2);
PhsData = AmpData;
for k = 1:cntrl(1)
    temp1 = textscan(fid,'%s %f %f',1);
    StnId{1,k} = temp1{1}{1};
    %StnId{1,k}  = T_o
    LatLon(k,:) = [temp1{2} temp1{3}];
    temp2 = textscan(fid,'%s %f %f %f %f',cntrl(2));
    AmpData{k,1}=temp2{2};
    AmpData{k,2}=temp2{3};
    PhsData{k,1}=temp2{4};
    PhsData{k,2}=temp2{5};
end;

fclose(fid);
Const = temp2{1};
NC = length(Const);
ststable = [];
errfid = fopen([runtitle,'_Outliers.txt'],'w');
stsname = {['R-Squared'],['Slope'],['Std-Dev'],['Avg-Err'],['Abs-Avg-Err'],['Norm-RMSE']};
stsvari = cell(0);
%% Get Color Palette and Specify Marker Size
mycolors = dlmread('harms.pal')/255;
%% Process Phase Constituents for Comparison
for k = 1:cntrl(1)
    obstest = PhsData{k,1};
    modtest = PhsData{k,2};
    for j = 1:cntrl(2)
        if obstest(j)~=flag && modtest(j)~=flag
            t1 = abs(obstest(j)-modtest(j));
            t2 = abs(-(360-obstest(j))-modtest(j));
            t3 = abs(obstest(j)+(360-modtest(j)));
            if min([t1 t2 t3]) == t2
                obstest(j)=-(360-obstest(j));
            elseif min([t1 t2 t3]) == t3
                modtest(j)=-(360-modtest(j));
            end
        end
    end
    PhsData{k,1}=obstest;
    PhsData{k,2}=modtest;
end

%% Process Bins
if isempty(bins) == 0
    bid = fopen(bins);
    temp = textscan(bid,'%d',1);
    NB = temp{1};
    for k = 1:NB
        temp = textscan(bid,'%s %d',1);
        mybins{k,1} = temp{1};
        temp1 = textscan(bid,'%s',temp{2});
        mybins{k,2} = temp1{1};
        for j = 1:temp{2}
            temp2 = find(strcmp(StnId,temp1{1}{j}));
            mybins{k,3}(j) = temp2;
        end
    end
end


%-------------------------------------------------------------------------%
%------------------------ Figure Creation --------------------------------%
%-------------------------------------------------------------------------%
%%  Plot Amplitudes and Phase 
ampobs = zeros(cntrl(1),1);
ampmod = ampobs;
phsobs = ampobs;
phsmod = ampobs;
allamp = [];
allphs = [];
for k = 1:NC
    for j = 1:cntrl(1)
        ampobs(j) = AmpData{j,1}(k);
        ampmod(j) = AmpData{j,2}(k);
        phsobs(j) = PhsData{j,1}(k);
        phsmod(j) = PhsData{j,2}(k);
        
        ampobs_org(j) = AmpData{j,1}(k);
        ampmod_org(j) = AmpData{j,2}(k);
        phsobs_org(j) = PhsData{j,1}(k);
        phsmod_org(j) = PhsData{j,2}(k);
    end
       
    rm = isnan(ampmod); 
    ampobs(rm)=[];
    ampmod(rm)=[];
    phsobs(rm)=[]; 
    phsmod(rm)=[];
    allamp = [allamp; ampobs ampmod];
    allphs = [allphs; phsobs phsmod];
    %LatLon(rm,:) =[];
   idx = find(~rm); 
   plotconst(ampobs,ampmod,[],k,[Const{k},' Amplitude ',runtitle],...
       [filetitle,'_',Const{k},'_amp'],[Const{k},'-Amplitude'],myMarks{k},1,0,ampaxis(k,:),'Amplitude')
   plotconst(phsobs,phsmod,20,NC+k,[Const{k},' Phase ',runtitle],...
       [filetitle,'_',Const{k},'_phs'],[Const{k},'-Phase'],myMarks{k},1,0,phsaxis(k,:),'Phase')
    if xml == 1
        createxmlfile(StnId,LatLon(:,:),ampobs_org,ampmod_org,...
            chop(nanstd(ampobs_org),1)/2,[Const{k},'_amp'],filetitle);
        createxmlfile(StnId,LatLon(:,:),phsobs_org,phsmod_org,20,...
            [Const{k},'_phs'],filetitle);
    end
end
%% Set-Up and Process Bins
indiv = plotbins;
if isempty(bins) == 0
    for k = 1:NB %length(bins(:,1))
        ampobs = [];
        ampmod = [];
        phsobs = [];
        phsmod = [];
        for j = 1:length(mybins{k,2})
           ampobs = [ampobs;AmpData{mybins{k,3}(j),1}(:)];
           ampmod = [ampmod;AmpData{mybins{k,3}(j),2}(:)];
           phsobs = [phsobs;PhsData{mybins{k,3}(j),1}(:)];
           phsmod = [phsmod;PhsData{mybins{k,3}(j),2}(:)];
        end
        markervec = cell(length(ampobs),1);
        for j = 1:NC
        markervec(j:NC:end) = myMarks(j);
        end
        plotconst(ampobs,ampmod,[],2*NC+k,[mybins{k,1}{1},' Amplitude '...
            ,runtitle],[filetitle,'_',mybins{k,1}{1},'_amp'],[mybins{k,1}{1},'-Amplitude'],markervec,0,2,ampext,'Amplitude')
        plotconst(phsobs,phsmod,20,3*NC+k,[mybins{k,1}{1},' Phase '...
            ,runtitle],[filetitle,'_',mybins{k,1}{1},'_phs'],[mybins{k,1}{1},'-Phase'],markervec,0,2,phsext,'Phase')
    end
end

%% Create Figures for all Amp and Phs data
indiv = plotall;
ampobs = allamp(:,1);
ampmod = allamp(:,2);
phsobs = allphs(:,1);
phsmod = allphs(:,2);
markervec = cell(length(ampobs),1);
cnt = 0;
for k = 1:NC
    for j = 1:cntrl(1)
        cnt = cnt+1;
        markervec(cnt) = myMarks(k);
    end
end

disp('creating all figs')

plotconst(ampobs(ampobs~=9999),ampmod(ampmod~=9999),[],4*NC+1,['All Amplitude ',runtitle],...
   [filetitle,'_','All_amp'],'All-Amplitude',markervec,0,2,ampext,'Amplitude')
plotconst(phsobs,phsmod,20,4*NC+2,['All Phase ',runtitle],...
   [filetitle,'_','All_phs'],'All-Phase',markervec,0,2,phsext,'Phase')

%% Construct Text File with desired Statistics
if sflag == 1;
fid = fopen([filetitle,'_stats_table.txt'],'w');
fprintf(fid,'%s\t','Constituent');
for k = 1:length(stsused)
    if stsused(k)==1
        fprintf(fid,'%s\t',stsname{k});
    end
end
fprintf(fid,'\n');
myloc = find(ststable(:,1) ~= 0);
for k = 1:length(myloc)
    fprintf(fid,'%s \t',stsvari{myloc(k)});
    for i = 1:length(ststable(1,:))
        fprintf(fid,'%f \t',ststable(myloc(k),i));
    end
    fprintf(fid,'\n');
end
fclose(fid);
end

fclose('all');
%-------------------------------------------------------------------------%
%----------- Internal Function for construction of Figures ---------------%
%-------------------------------------------------------------------------%
function plotconst(obs,mod,errband,fignum,figtitle,filename,stsname,mark,plotopt,errrec,myaxis,dtype)
global mycolors flag ms ststable stsused filetypes stsvari errfid StnId indiv fsz fszstat
% Remove any null-data points (if desired)
loc = find(obs~=flag & mod~=flag);
obs = obs(loc);
mod = mod(loc);
if plotopt ~= 1
    mark = mark(loc);
end
if errrec == 0
mylocstn = StnId(loc);
end
if isempty(errband)
errband = chop(std(obs),1)/2;
end
% Determine Data range for Figure
mymin = floor(min([min(obs) min(mod)])*1000)/1000;
mymax = ceil(max([max(obs) max(mod)])*1000)/1000;
mymin = min(mymin,myaxis(1));
mymax = max(mymax,myaxis(2));
%-------------------- Statistical Analysis -------------------------------%
% R-Squared Value
%R2 = (sum((obs-mean(obs)).*(mod-mean(mod)))/(sqrt(sum((obs-mean(obs)).^2))*sqrt(sum((mod-mean(mod)).^2))))^2;
% R2 = (corr(obs,mod))^2;
obs(isnan(obs))=[]; mod(isnan(mod))=[];
R2 = (sum((obs-nanmean(obs)).*(mod-nanmean(mod)))/(sqrt(sum((obs-nanmean(obs)).^2))*sqrt(sum((mod-nanmean(mod)).^2))))^2;
% SStot = sum((mod-mean(mod)).^2);
% SSreg = sum((obs-mean(mod)).^2);
% SSerr = sum((mod-obs).^2);
% R2 = 1-SSerr/SStot;
% Best Fitting Line
options = optimset('Display','off','TolFun',1e-10);
a= lsqcurvefit(@(a,obs) a*obs,1,obs,mod,[],[],options);
%
MOerr = mod-obs;
% Standard Deviation
mystd = std(MOerr);
% Average Difference between Observed and Modeled
avgerr = mean(MOerr);
% Average of the Absolute Difference
avgabs = mean(abs(MOerr));
% Normalized RMSE
nrmse = sqrt(sum((mod-obs).^2)./sum(obs.^2));
% Adendum Stats Table
sts = [R2 a mystd avgerr avgabs nrmse];
ststable(fignum,:) = sts(stsused==1);
stsvari(fignum) = {stsname};
%------------------------ Error Bands ------------------------------------%
% Construct Error Bands
e0 = obs;
u1 = e0+errband;
u2 = e0+2*errband;
u3 = e0+3*errband;
l1 = e0-errband;
l2 = e0-2*errband;
l3 = e0-3*errband;
% Determine data within each errorband
du0 = find(((mod-e0)>=0).*((mod-u1)<0));
du1 = find(((mod-u1)>=0).*((mod-u2)<0));
du2 = find(((mod-u2)>=0).*((mod-u3)<0));
du3 = find(((mod-u3)>=0));
dl0 = find(((e0-mod)>=0).*((l1-mod)<0));
dl1 = find(((l1-mod)>=0).*((l2-mod)<0));
dl2 = find(((l2-mod)>=0).*((l3-mod)<0));
dl3 = find((l3-mod)>=0);
% Recreate Error Bands for figure
e0 = [mymin mymax];
u1 = e0+errband;
u2 = e0+2*errband;
u3 = e0+3*errband;
l1 = e0-errband;
l2 = e0-2*errband;
l3 = e0-3*errband;

% Record major outliers
if isempty(du3)+isempty(dl3)+errrec<2
    fprintf(errfid,'%s \n',figtitle);
    for i = 1:length(du3)
        fprintf(errfid,'%s \n',mylocstn{du3(i)});
        %fprintf(errfid,'%s \n',T_o.Name{du3(i)});
        
    end
    for i = 1:length(dl3)
        fprintf(errfid,'%s \n',mylocstn{dl3(i)});
       %fprintf(errfid,'%s \n',T_o.Name{dl3(i)});
    end
end

%-------------------- Construct Figure -----------------------------------%
if indiv == 1
figure(fignum)
clf; set(gcf,'renderer','painters');set(gcf,'color','white'); figure(fignum);        
    set(gcf,'PaperType','usLetter','PaperOrientation','portrait',...
            'PaperUnits','Inches','PaperPositionMode','Manual') 
    set(gcf,'PaperPosition',[0 0 3.5 3.5],'Position',[100 100 350 350])
if plotopt == 1
    hold on
    plot(obs(du0),mod(du0),'Marker',mark,'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(5,:),'MarkerFaceColor',mycolors(5,:))
    plot(obs(du1),mod(du1),'Marker',mark,'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(6,:),'MarkerFaceColor',mycolors(6,:))
    plot(obs(du2),mod(du2),'Marker',mark,'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(7,:),'MarkerFaceColor',mycolors(7,:))
    plot(obs(du3),mod(du3),'Marker',mark,'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(8,:),'MarkerFaceColor',mycolors(8,:))
    plot(obs(dl0),mod(dl0),'Marker',mark,'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(4,:),'MarkerFaceColor',mycolors(4,:))
    plot(obs(dl1),mod(dl1),'Marker',mark,'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(3,:),'MarkerFaceColor',mycolors(3,:))
    plot(obs(dl2),mod(dl2),'Marker',mark,'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(2,:),'MarkerFaceColor',mycolors(2,:))
    plot(obs(dl3),mod(dl3),'Marker',mark,'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(1,:),'MarkerFaceColor',mycolors(1,:))
else
    hold on
    for n = 1:length(du0)
        plot(obs(du0(n)),mod(du0(n)),'Marker',mark{du0(n)},'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(5,:),'MarkerFaceColor',mycolors(5,:))
    end
    for n = 1:length(du1)
        plot(obs(du1(n)),mod(du1(n)),'Marker',mark{du1(n)},'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(6,:),'MarkerFaceColor',mycolors(6,:))
    end
    for n = 1:length(du2)
        plot(obs(du2(n)),mod(du2(n)),'Marker',mark{du2(n)},'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(7,:),'MarkerFaceColor',mycolors(7,:))
    end
    for n = 1:length(du3)
        plot(obs(du3(n)),mod(du3(n)),'Marker',mark{du3(n)},'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(8,:),'MarkerFaceColor',mycolors(8,:))
    end
    for n = 1:length(dl0)
        plot(obs(dl0(n)),mod(dl0(n)),'Marker',mark{dl0(n)},'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(4,:),'MarkerFaceColor',mycolors(4,:))
    end
    for n = 1:length(dl1)
        plot(obs(dl1(n)),mod(dl1(n)),'Marker',mark{dl1(n)},'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(3,:),'MarkerFaceColor',mycolors(3,:))
    end
    for n = 1:length(dl2)
        plot(obs(dl2(n)),mod(dl2(n)),'Marker',mark{dl2(n)},'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(2,:),'MarkerFaceColor',mycolors(2,:))
    end
    for n = 1:length(dl3)
        plot(obs(dl3(n)),mod(dl3(n)),'Marker',mark{dl3(n)},'LineStyle','none','MarkerSize',ms,'MarkerEdgeColor',mycolors(1,:),'MarkerFaceColor',mycolors(1,:))
    end
end
if stsused(2) == 1
    plot([mymin mymax],a*[mymin mymax],'linewidth',2)
end
plot(e0,e0,'--',e0,u1,e0,u2,e0,u3,e0,l1,e0,l2,e0,l3,'Color',mycolors(9,:))
hold off
title(figtitle,'fontsize',fsz,'interpreter','latex')
xlabel(['Observed ',dtype],'fontsize',fsz,'interpreter','latex')
ylabel(['Model ',dtype],'fontsize',fsz,'interpreter','latex')
axis([mymin mymax mymin mymax]);
% axis square
%--------------- Include Statistics on Figure ----------------------------%
stslab = {['\(R^2=\)',num2str(R2)],['\(y=\)',num2str(a),'x'],['\(\sigma=\)',num2str(mystd)],['\(\bar{\epsilon}=\)',num2str(avgerr)],['\(\bar{\left|\epsilon\right|}=\)',num2str(avgabs)],['\(E=\)',num2str(nrmse)]};
stslab = stslab(stsused==1);
stslab = [stslab ['\(\Delta=\)',num2str(errband)]];
rng = mymax-mymin;
text(rng/10*7.5+mymin,rng/10*2+mymin,stslab,'interpreter','latex','fontsize',fszstat);
for k = 1:length(filetypes)
%     saveas(gcf,[filename,filetypes{k,1}],filetypes{k,2});
    print(filetypes{k},'-r300',filename);
end
end
%-------------------------------------------------------------------------%
%----------- Internal Function for construction of XML Files -------------%
%-------------------------------------------------------------------------%
function createxmlfile(StnId,LatLon,obs,mod,delta,xmlabel,runtitle)
global flag
obs(obs==flag) = NaN;
mod(mod==flag) = NaN;
fid = fopen([runtitle,'_',xmlabel,'.xml'],'w');
fprintf(fid,'%s \n','<markers>');
for k = 1:length(LatLon)
    str = ['<marker stnid="',StnId{k},'" lat="',num2str(LatLon(k,2)),...
        '" lng="',num2str(LatLon(k,1)),'" label="',xmlabel,...
        '" observed="',num2str(obs(k)),'" model="',num2str(mod(k)),'" error="',num2str(-obs(k)+mod(k)),...
        '" delta="',num2str(delta),'"/>'];
    fprintf(fid,'%s \n',str);
end
fprintf(fid,'%s','</markers>');
fclose(fid);
