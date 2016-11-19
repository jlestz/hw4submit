% function for plotting final Temperature 
% data created by heat_serial, heat_omp, heat_mpi 
% opens files of form heat_<str>.<nx>.<np>.all
% (.<np>.all omitted for serial results)

function [T,Tavg,trun] = plot_temp(str,nx,np,do_save)

% by default, do not save data and close the figure
if nargin < 4
    do_save=0; 
end

% construct file name and title for plot 
fname=['out_' str '.' num2str(nx)];
tstr=[str ' nx = ' num2str(nx)];
if ~strcmp(str,'serial')
    fname=[fname '.' num2str(np) '.all'];
    tstr=[tstr ' nprocs = ' num2str(np)];
end
tstr=[tstr ' T_{avg} = '];

% open the file 
f=fopen(fname);

% store the final temperature data 
T=fscanf(f,'%f',[nx nx]);

% read extraneous information 
fgetl(f); 
fgetl(f); 

% get the run time 
stime=fgetl(f); 
split=strsplit(stime,'='); 
trun=str2double(split{2}); 

% get the average temperature 
s=fgetl(f);
sdex=strfind(s,'0');
s=s(sdex:end);
Tavg=str2double(s); 
tstr=[tstr s];

% close the file 
fclose(f);

% set the grid for plotting 
xpts=linspace(0,pi,nx+1);
xpts=xpts(1:end-1);
ypts=linspace(0,pi,nx+2);
ypts=ypts(2:end-1);

% normalize 
xpts=xpts/pi; 
ypts=ypts/pi; 

% create a 2D mesh
[X,Y]=meshgrid(xpts,ypts);

% plot the temperature data 
FS=14;
figure ;
h=pcolor(X,Y,T); 
set(h,'linestyle','none');
colorbar;
caxis([0 1]); 
set(gca,'fontsize',FS);
box on ;
title(tstr);
axis tight

% save the figure to file 
if do_save 
    hgexport(gcf, [fname '.png'], hgexport('factorystyle'), 'Format', 'png');
    close(gcf); 
end

end
