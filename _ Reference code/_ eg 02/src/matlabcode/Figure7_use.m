clear all
close all
clc


h=figure(1);
set(h, 'position', [50   50   800   680]);

bb=12;
cc = [-2:1:50];

%% 
han01 = subplot(3,2,1);

load data/plot2_data_1;

m_proj('Mollweide','lon',[20 380], 'lat',[-90 90]);
m_grid('xtick',[20:360:380], 'xticklabel', '', 'ytick',[-90:30:90], 'yticklabel', '', 'tickdir','out','FontName', 'Times New Roman', 'fontsize', 10, 'linestyle', ':');
hold on;
[C1, h1] = m_contourf (X, Y, Z2, cc, 'LineStyle', 'none' );
m_gshhs_i('patch',[.7 0.7 .7]);

caxis([0 bb]);
clear('X', 'Y', 'Z1', 'Z2');

%% 
han02 = subplot(3,2,2);

load data/plot2_data_2;

m_proj('Mollweide','lon',[20 380], 'lat',[-90 90]);
m_grid('xtick',[20:360:380], 'xticklabel', '', 'ytick',[-90:30:90],'yticklabel', '','tickdir','out','FontName', 'Times New Roman', 'fontsize', 10, 'linestyle', ':');
hold on;
[C1, h1] = m_contourf (X, Y, Z2, cc, 'LineStyle', 'none' );
m_gshhs_i('patch',[.7 0.7 .7]);

caxis([0 bb]);
clear('X', 'Y', 'Z1', 'Z2');


%% han04, 2000
han03 = subplot(3,2,3);

load data/plot2_data_3;

m_proj('Mollweide','lon',[20 380], 'lat',[-90 90]);
m_grid('xtick',[20:360:380], 'xticklabel', '', 'ytick',[-90:30:90],'yticklabel', '','tickdir','out','FontName', 'Times New Roman', 'fontsize', 10, 'linestyle', ':');
hold on;
[C1, h1] = m_contourf (X, Y, Z2, cc, 'LineStyle', 'none' );
m_gshhs_i('patch',[.7 0.7 .7]);

caxis([0 bb]);
clear('X', 'Y', 'Z1', 'Z2');

%% 
han04 = subplot(3,2,4);

load data/plot2_data_4;

m_proj('Mollweide','lon',[20 380], 'lat',[-90 90]);
m_grid('xtick',[20:360:380], 'xticklabel', '', 'ytick',[-90:30:90],'yticklabel', '','tickdir','out','FontName', 'Times New Roman', 'fontsize', 10, 'linestyle', ':');
hold on;
[C1, h1] = m_contourf (X, Y, Z2, cc, 'LineStyle', 'none' );
m_gshhs_i('patch',[.7 0.7 .7]);

caxis([0 bb]);
clear('X', 'Y', 'Z1', 'Z2');



%% 

han05=subplot(3,2,5);

load data/plot2_data_5;

m_proj('Mollweide','lon',[20 380], 'lat',[-90 90]);
m_grid('xtick',[20:360:380], 'xticklabel', '', 'ytick',[-90:30:90],'yticklabel', '','tickdir','out','FontName', 'Times New Roman', 'fontsize', 10, 'linestyle', ':');
hold on;
[C1, h1] = m_contourf (X, Y, Z2, cc, 'LineStyle', 'none' );
m_gshhs_i('patch',[.7 0.7 .7]);

caxis([0 bb]);
clear('X', 'Y', 'Z1', 'Z2');

%% 
han06 = subplot(3,2,6);

load data/plot2_data_6;

m_proj('Mollweide','lon',[20 380], 'lat',[-90 90]);
m_grid('xtick',[20:360:380], 'xticklabel', '', 'ytick',[-90:30:90],'yticklabel', '','tickdir','out','FontName', 'Times New Roman', 'fontsize', 10, 'linestyle', ':');
hold on;
[C1, h1] = m_contourf (X, Y, Z2, cc, 'LineStyle', 'none' );
m_gshhs_i('patch',[.7 0.7 .7]);

caxis([0 bb]);
clear('X', 'Y', 'Z1', 'Z2');

%%
colormap(parula());

han_cb=colorbar('southoutside', 'YTick', [0:2:bb], 'FontName', 'Times New Roman', 'FontSize', 14, 'TickDirection', 'out');
% set(han_cb, 'YTickLabel', {'7.70' '7.75' '7.80' '7.85' '7.90' '7.95' '8.00' '8.05' '8.10' '8.15' '8.20' '8.25'});
set(han_cb, 'position', [0.16    0.12    0.7    0.035]);

% text(4.4, 3.04, 'pH on Total Scale',  'FontName', 'Times New Roman', 'FontSize', 14, 'Rotation', 270);


%% setting positions:
A = get(han01,'position');
R = 1.1;
set(han01, 'position', [0.12   0.70  A(3)*R  A(4)*R]);
set(han02, 'position', [0.54   0.70  A(3)*R  A(4)*R]);
set(han03, 'position', [0.12   0.43  A(3)*R  A(4)*R]);
set(han04, 'position', [0.54   0.43  A(3)*R  A(4)*R]);
set(han05, 'position', [0.12   0.16  A(3)*R  A(4)*R]);
set(han06, 'position', [0.54   0.16  A(3)*R  A(4)*R]);

% text:
text(-4.20, 6.20,  'A) 1800 climates disappearing in 2000', 'FontName', 'Times New Roman', 'Fontsize', 14); 
text( 0.3,  6.20,  'B) Novel climates in 2000 since 1800', 'FontName', 'Times New Roman', 'Fontsize', 14); 
text(-4.20, 3.70,  'C) 2000 climates disappearing in 2100 RCP 4.5', 'FontName', 'Times New Roman', 'Fontsize', 14); 
text( 0.3,  3.70,  'D) Novel climates in 2100 RCP 4.5 since 2000', 'FontName', 'Times New Roman', 'Fontsize', 14); 
text(-4.20, 1.20,  'E) 2000 climates disappearing in 2100 RCP 8.5', 'FontName', 'Times New Roman', 'Fontsize', 14); 
text( 0.3,  1.20,  'F) Novel climates in 2100 RCP 8.5 since 2000', 'FontName', 'Times New Roman', 'Fontsize', 14); 


text( -1, -2.0,  'Mahalanobis distance', 'FontName', 'Times New Roman', 'Fontsize', 16); 











