function syncAxes(handles)

% synchronize axes
linkaxes([handles.axes_surf, handles.axes_spectrum], 'x');
addlistener( handles.axes_surf, 'YLim', 'PostSet', @(src,evt) y2xlimsync(src, evt, handles.axes_kinetic, handles.axes_surf));
addlistener( handles.axes_kinetic, 'XLim', 'PostSet', @(src,evt) x2ylimsync(src, evt, handles.axes_surf, handles.axes_kinetic));

end

