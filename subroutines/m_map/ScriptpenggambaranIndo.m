%Load Garis Pantai
Garpan=load('GarpanPakIbnu.dat');
Longi=Garpan(:,1);
Latti=Garpan(:,2);
%% Penggambaran
%plot
figure(1);
m_proj('miller','long',[40 120],'lat',[-15 15]);
m_coast;
m_plot(Longi,Latti,'black','LineWidth',1);hold on;
m_patch(Longi,Latti,[0.4 0.4 0.4]);
m_grid('box','fancy','tickdir','in','xtick',7,'ytick',15,'linestyle','none','fontsize',15,'linewidth',3)
%Sumbu-x,Sumbu-y,SkalaWarna
%ax=gca;
%set(ax,'FontSize',15,'FontWeight','bold');
%cbh=colorbar('FontWeight','bold','FontName','Times News Roman');
%Colorbar title
xlabel('Longitude','FontWeight','bold','FontName','Times News Roman','Fontsize',15);
ylabel('Latitude','FontWeight','bold','FontName','Times News Roman','Fontsize',15);
title(strcat('Peta Perairan Indonesia'),'Fontweight','bold','Fontsize',15);
%SavingFigureFullScreen
%screen_size=get(0,'ScreenSize');
%origSize=get(figure(1),'Position');%grab original on screen size
%set(figure(1),'Position',[0 0 screen_size(3) screen_size(4) ] );%set tp screen size 
%set(figure(1),'PaperPositionMode','auto');%set paper pos for printing
%saveas(figure(1),strcat('Average Wave Height25years.png'));
%set(figure(1),'Position',origSize);%set back to original dimensions
%close all;
