function refresh1DPlots(handles, wCursor, tIndCursor)

TA = handles.TA;

tIndCursor = round(tIndCursor);
[~, wInd] = min(abs(TA.WVec - wCursor));

% draw point under cursor
% % % if (handles.checkbox_kinetic_semilogx.Value == 1)
% % %     semilogx(handles.axes_kinetic, TA.TVec(TA.TVec>0), TA.deltaA(TA.TVec>0, wInd)*1e3, 'Color', 'k', 'LineWidth', 0.5);
% % % else
% % %     plot(handles.axes_kinetic, TA.TVec, TA.deltaA(:, wInd)*1e3, 'Color', 'k', 'LineWidth', 0.5);
% % % end
plot(handles.axes_kinetic_lin, TA.TVec(TA.TVec<2), TA.deltaA(TA.TVec<2, wInd)*1e3, 'Color', 'k', 'LineWidth', 0.5);
semilogx(handles.axes_kinetic_log, TA.TVec(TA.TVec>2), TA.deltaA(TA.TVec>2, wInd)*1e3, 'Color', 'k', 'LineWidth', 0.5);
plot(handles.axes_spectrum, TA.WVec, TA.deltaA(tIndCursor, :)*1e3, 'Color', 'k', 'LineWidth', 0.5);

% draw saved points-of-interest
hold(handles.axes_kinetic_lin, 'on')
hold(handles.axes_kinetic_log, 'on')
hold(handles.axes_spectrum, 'on')
poiList = handles.poiList;
for i=1:poiList.Num_nodes
    p = poiList.get_node(i);
    [~, w] = min(abs(p.w-TA.WVec));
% % %     if (handles.checkbox_kinetic_semilogx.Value == 1)
% % %         semilogx(handles.axes_kinetic, TA.TVec(TA.TVec>0), TA.deltaA(TA.TVec>0, w)*1e3, 'Color', p.Color, 'LineWidth', 1);
% % %     else
% % %         plot(handles.axes_kinetic, TA.TVec, TA.deltaA(:, w)*1e3, 'Color', p.Color, 'LineWidth', 1);
% % %     end
    plot(handles.axes_kinetic_lin, TA.TVec(TA.TVec<2), TA.deltaA(TA.TVec<2, w)*1e3, 'Color', p.Color, 'LineWidth', 1);
    semilogx(handles.axes_kinetic_log, TA.TVec(TA.TVec>2), TA.deltaA(TA.TVec>2, w)*1e3, 'Color', p.Color, 'LineWidth', 1);
    plot(handles.axes_spectrum, TA.WVec, TA.deltaA(round(p.t), :)*1e3, 'Color', p.Color, 'LineWidth', 1);
end
hold(handles.axes_kinetic_lin, 'off')
hold(handles.axes_kinetic_log, 'off')
hold(handles.axes_spectrum, 'off')

% set axes lalbels
if handles.fsTA_radio.Value == 1
    xlabel(handles.axes_kinetic_log, 'Time (ps)');
else
    xlabel(handles.axes_kinetic_log, 'Time (\mus)');
end

ylabel(handles.axes_kinetic_lin, '\DeltaA (mOD)');
xlabel(handles.axes_spectrum, 'Wavelength (nm)');
ylabel(handles.axes_spectrum, '\DeltaA (mOD)');

% set axes limits
% % % if (handles.checkbox_kinetic_semilogx.Value == 1)
% % %     xlim(handles.axes_kinetic, [0.1, max(TA.TVec(TA.TVec>0))]);
% % % else
% % %     xlim(handles.axes_kinetic, [min(TA.TVec), max(TA.TVec)]);
% % % end
xlim(handles.axes_kinetic_lin, [min(TA.TVec), 2]);
xlim(handles.axes_kinetic_log, [2, max(TA.TVec)]);
ylim(handles.axes_kinetic_lin, [min(min(TA.deltaA)), max(max(TA.deltaA))]*1e3);
ylim(handles.axes_kinetic_log, [min(min(TA.deltaA)), max(max(TA.deltaA))]*1e3);
xlim(handles.axes_spectrum, [min(TA.WVec), max(TA.WVec)]);
ylim(handles.axes_spectrum, [min(min(TA.deltaA)), max(max(TA.deltaA))]*1e3);

% add legends
legend(handles.axes_kinetic_log, strcat(num2str(TA.WVec(wInd),'%.2f'), ' nm'), 'Location', 'northeast');
if handles.fsTA_radio.Value == 1
    legend(handles.axes_spectrum, strcat(num2str(TA.TVec(tIndCursor),'%.1f'), ' ps'), 'Location', 'northeast');
else
    legend(handles.axes_spectrum, strcat(num2str(TA.TVec(tIndCursor)*1e3,'%.1f'), ' ns'), 'Location', 'northeast');
end

% disable y-axis ticks on logarithmic kinetics plot
handles.axes_kinetic_log.YTickLabel = '';

end

