% function to plot all temperature data, looping over plot_temp.m
% plots

function [Tavgall,trunall] = plot_temp_all(modes,nx,np,do_save,do_summary)

if nargin < 2 
    do_summary=0; 
end 

if nargin < 1 
    do_save=0; 
end

% exhaustive list of run parameters
% modes={'serial','omp','mpi'};
% nx=[128 256 512];
% np=[1 2 4 8 16];

% number of each parameter
nummodes=numel(modes);
numx=numel(nx);
nump=numel(np);

% for storing the average temperature and run times
Tavgall = zeros(nummodes,numx,nump);
trunall = Tavgall;

% loop over all parameters
for im=1:nummodes
    for ix=1:numx
        for ip=1:nump
            [~,Tavg,trun]=plot_temp(modes{im},nx(ix),np(ip),do_save);
            Tavgall(im,ix,ip)=Tavg;
            trunall(im,ix,ip)=trun;
        end
    end
end

if do_summary
    % exclude fakes from plotting
    Tavgall(1,:,2:end)=NaN;
    Tavgall(2,:,end)=NaN;
    
    % make plots of runtime dependency on np
    MS=10; FS=14;
    figure ;
    for j=1:numx
        subplot(1,numx,j);
        for i=1:nummodes
            semilogy(np,reshape(trunall(i,j,:),size(np)),'-ok','markersize',MS); hold on ;
        end
        set(gca,'fontsize',FS);
        legend('serial','omp','mpi');
        xlabel('Grid Size');
        ylabel('Run Time');
        title(['nx = ' num2str(nx(j))]);
    end
    
    % make plots of runtime dependencey on nx
    figure ;
    for j=1:nump
        subplot(1,nump,j);
        for i=1:nummodes
            semilogy(nx,reshape(trunall(i,:,j),size(nx)),'-ok','markersize',MS); hold on ;
        end
        set(gca,'fontsize',FS);
        legend('serial','omp','mpi');
        xlabel('Number of Processors');
        ylabel('Run Time');
        title(['np = ' num2str(np(j))]);
    end
end


end