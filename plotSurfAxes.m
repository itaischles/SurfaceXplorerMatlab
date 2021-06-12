function plotSurfAxes(handles)

TA = handles.TA;

% scale TA surface
% % % scaled_TA = (TA.deltaA - min(min(TA.deltaA)));
% % % scaled_TA = scaled_TA./max(max(scaled_TA));
scaled_TA = atan(TA.deltaA*1000);

% plot scaled TA surface
% % % climits = [str2num(handles.lowlim_text.String), str2num(handles.uplim_text.String)];
% % % levellist = linspace(climits(1),climits(2),20);
% % % imagesc(handles.axes_surf, [min(TA.WVec),max(TA.WVec)], [1, numel(TA.TVec)], TA.deltaA, climits);
% % % contourf(handles.axes_surf, TA.WVec, 1:1:numel(TA.TVec), scaled_deltaA, 'LineColor', 'none', 'LevelList', levellist);
TAsurf = pcolor(handles.axes_surf, TA.WVec, 1:1:numel(TA.TVec), scaled_TA);
set(TAsurf, 'EdgeColor', 'none');
shading interp

% set color map
% % % cmap = bluewhitered(21, handles.axes_surf.CLim);
cmap = flipud(cbrewer('RdBu', 128, 'pchip'));
% % % cmap = jet(20);
colormap(handles.axes_surf, cmap)

% % % colormap(handles.axes_surf, jet)



% set axes properties
handles.axes_surf.YTick = [];
% % % handles.axes_surf.XLim = [min(TA.WVec),max(TA.WVec)];
handles.axes_surf.YDir = 'normal';

% set axes labels
t = handles.TA.TVec;
tlabels = {'-2','0','10','100','1000'};
tinds = zeros(numel(tlabels),1);
for i=1:numel(tinds)
    [~,tinds(i)] = min(abs(t-str2num(tlabels{i})));
end
xlabel(handles.axes_surf, 'Wavelength (nm)');
ylabel(handles.axes_surf, 'Time (ps)');
handles.axes_surf.YTick = tinds;
handles.axes_surf.YTickLabel = tlabels;
handles.axes_surf.TickDir = 'out';

% draw points-of-interest
hold(handles.axes_surf, 'on')
poiList = handles.poiList;
for i=1:poiList.Num_nodes
    p = poiList.get_node(i);
    plot(handles.axes_surf, p.w, p.t, 'LineStyle', 'none', 'Marker', 'o', 'Color', p.Color, 'MarkerFaceColor', p.Color, 'MarkerSize', 4);
end
hold(handles.axes_surf, 'off')

% save
guidata(handles.fig_surfxplorer, handles);

end

