function updateGui(handles)

% reset points-of-interest array
handles.POI = cell(0,1);

% initialize axes
cla(handles.axes_surf, 'reset');
cla(handles.axes_kinetic_lin, 'reset');
cla(handles.axes_kinetic_log, 'reset');
cla(handles.axes_spectrum, 'reset');

% disable clear points button
handles.button_clearPoints.Enable = 'off';

% initialize lower and upper limits of colormap
% % % lowlim = min(min(handles.TA.deltaA));
% % % uplim = max(max(handles.TA.deltaA));
% % % handles.lowlim_text.String = num2str(lowlim*1000);
% % % handles.uplim_text.String = num2str(uplim*1000);
handles.lowlim_text.String = num2str(0);
handles.uplim_text.String = num2str(1);

% save
guidata(handles.fig_surfxplorer, handles);

end

