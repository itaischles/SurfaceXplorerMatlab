function plotFit(handles, fit, fs_or_ns, model_name)

mainfig = figure('Position', [192 72 980 591]);

if strcmp(fs_or_ns, 'fsTA')
    ind = handles.TA.TVec < 1;
else
    ind = handles.TA.TVec < 0.01;
end
ax1 = axes(mainfig, 'Position', [0.5552 0.0788 0.0695 0.3901], 'Units', 'normalized');
plot(handles.TA.TVec, fit.Population/max(max(fit.Population)), 'LineWidth', 2);
ylabel('Normalized model population')
xlim([min(handles.TA.TVec) max(handles.TA.TVec(ind))])
ax2 = axes(mainfig, 'Position', [0.6247 0.0788 0.3475 0.3901], 'Units', 'normalized');
plot(handles.TA.TVec, fit.Population/max(max(fit.Population)), 'LineWidth', 2)
xlim([min(handles.TA.TVec(~ind)), max(handles.TA.TVec(~ind))])
set(gca,'Xscale','log','YTickLabel','')
linkaxes([ax1 ax2],'y')
if strcmp(fs_or_ns, 'fsTA')
    xlabel('Time (ps)')
else
    xlabel('Time (\mus)')
end
switch model_name
    case {'A->B', 'A->B->GND', '2A->B', '2A->B->GND'}
        legend('A','B', 'Location', 'best');
    case {'A->B->C', 'A->B->C->GND', 'A->B,2B->C'}
        legend('A','B','C', 'Location', 'best');
    case {'A->B->C->D'}
        legend('A','B','C','D', 'Location', 'best');
end
legend('boxoff');
set([ax1,ax2],...
    'Box', 'on',...
    'LineWidth', 1.2,...
    'XMinorTick', 'on',...
    'YMinorTick', 'on',...
    'TickLength', [0.007 0.007],...
    'FontName', 'Helvetica',...
    'FontSize', 10);

if strcmp(fs_or_ns, 'fsTA')
    time_unit = ' ps';
else
    time_unit = '\mus';
end
ax3 = axes(mainfig, 'Position', [0.5552 0.5656 0.4170 0.3901], 'Units', 'normalized');
plot(ax3, handles.TA.WVec, zeros(numel(handles.TA.WVec), 1), 'Color', [0.2,0.2,0.2]);
ax3.ColorOrderIndex = 1;
hold on
h = plot(handles.TA.WVec, fit.SpectraFit*1000, 'LineWidth', 2);
xlabel('Wavelength (nm)')
ylabel('\DeltaA (mOD)')
xlim([min(handles.TA.WVec), max(handles.TA.WVec)])
switch model_name
    case 'A->B'
        tauAB_txt = strcat('    \tau_{A->B}=', num2str(fit.FitParams(end), '%.2f'), '\pm', num2str(fit.Std(end), '%.2f'), time_unit);
        legend(h, strcat('A ',tauAB_txt), 'B', 'Location', 'best');
    case 'A->B->C'
        tauAB_txt = strcat('    \tau_{A->B}=', num2str(fit.FitParams(end-1), '%.2f'), '\pm', num2str(fit.Std(end-1), '%.2f'), time_unit);
        tauBC_txt = strcat('    \tau_{B->C}=', num2str(fit.FitParams(end), '%.2f'), '\pm', num2str(fit.Std(end), '%.2f'), time_unit);
        legend(h, strcat('A ',tauAB_txt), strcat('B ',tauBC_txt), 'C', 'Location', 'best');
    case 'A->B->C->D'
        tauAB_txt = strcat('    \tau_{A->B}=', num2str(fit.FitParams(end-2), '%.2f'), '\pm', num2str(fit.Std(end-2), '%.2f'), time_unit);
        tauBC_txt = strcat('    \tau_{B->C}=', num2str(fit.FitParams(end-1), '%.2f'), '\pm', num2str(fit.Std(end-1), '%.2f'), time_unit);
        tauCD_txt = strcat('    \tau_{C->D}=', num2str(fit.FitParams(end), '%.2f'), '\pm', num2str(fit.Std(end), '%.2f'), time_unit);
        legend(h, strcat('A ',tauAB_txt), strcat('B ',tauBC_txt), strcat('C ',tauCD_txt), 'D', 'Location', 'best');
    case 'A->B->GND'
        tauAB_txt = strcat('    \tau_{A->B}=', num2str(fit.FitParams(end-1), '%.2f'), '\pm', num2str(fit.Std(end-1), '%.2f'), time_unit);
        tauBG_txt = strcat('    \tau_{B->GND}=', num2str(fit.FitParams(end), '%.2f'), '\pm', num2str(fit.Std(end), '%.2f'), time_unit);
        legend(h, strcat('A ',tauAB_txt), strcat('B ',tauBG_txt), 'Location', 'best');
    case 'A->B->C->GND'
        tauAB_txt = strcat('    \tau_{A->B}=', num2str(fit.FitParams(end-2), '%.2f'), '\pm', num2str(fit.Std(end-2), '%.2f'), time_unit);
        tauBC_txt = strcat('    \tau_{B->C}=', num2str(fit.FitParams(end-1), '%.2f'), '\pm', num2str(fit.Std(end-1), '%.2f'), time_unit);
        tauCG_txt = strcat('    \tau_{C->GND}=', num2str(fit.FitParams(end), '%.2f'), '\pm', num2str(fit.Std(end), '%.2f'), time_unit);
        legend(h, strcat('A ',tauAB_txt), strcat('B ',tauBC_txt), strcat('C ',tauCG_txt), 'Location', 'best');
    case '2A->B'
        tau2AB_txt = strcat('    \tau_{2A->B}=', num2str(fit.FitParams(end), '%.2f'), '\pm', num2str(fit.Std(end), '%.2f'), time_unit);
        legend(h, strcat('A ',tau2AB_txt), 'B', 'Location', 'best');
    case 'A->B,2B->C'
        tauAB_txt = strcat('    \tau_{A->B}=', num2str(fit.FitParams(end-1), '%.2f'), '\pm', num2str(fit.Std(end-1), '%.2f'), time_unit);
        tau2BC_txt = strcat('    \tau_{2B->C}=', num2str(fit.FitParams(end), '%.2f'), '\pm', num2str(fit.Std(end), '%.2f'), time_unit);
        legend(h, strcat('A ',tauAB_txt), strcat('B ',tau2BC_txt), 'C', 'Location', 'best');     
    case '2A->B->GND'
        tau2AB_txt = strcat('    \tau_{2A->B}=', num2str(fit.FitParams(end-1), '%.2f'), '\pm', num2str(fit.Std(end-1), '%.2f'), time_unit);
        tauBG_txt = strcat('    \tau_{B->G}=', num2str(fit.FitParams(end), '%.2f'), '\pm', num2str(fit.Std(end), '%.2f'), time_unit);
        legend(h, strcat('A ',tau2AB_txt), strcat('B ',tauBG_txt), 'Location', 'best');     
end
legend('boxoff');
set(gca,...
    'Box', 'on',...
    'LineWidth', 1.2,...
    'XMinorTick', 'on',...
    'YMinorTick', 'on',...
    'TickLength', [0.007 0.007],...
    'FontName', 'Helvetica',...
    'FontSize', 10);

if strcmp(fs_or_ns, 'fsTA')
    ind = handles.TA.TVec < 1;
else
    ind = handles.TA.TVec < 0.01;
end
ax4 = axes(mainfig, 'Position', [0.0513 0.0788 0.0695 0.3901], 'Units', 'normalized');
plot(ax4, handles.TA.TVec, zeros(numel(handles.TA.TVec), 1), 'Color', [0.2,0.2,0.2]);
ax4.ColorOrderIndex = 1;
hold on
xlim([min(handles.TA.TVec) max(handles.TA.TVec(ind))]);
ax5 = axes(mainfig, 'Position', [0.1208 0.0788 0.3649 0.3901], 'Units', 'normalized');
plot(ax5, handles.TA.TVec, zeros(numel(handles.TA.TVec), 1), 'Color', [0.2,0.2,0.2]);
ax5.ColorOrderIndex = 1;
hold on
xlim([min(handles.TA.TVec(~ind)), max(handles.TA.TVec(~ind))])
set(gca,'Xscale','log','YTickLabel','')
linkaxes([ax4 ax5],'y')
h = zeros(1, numel(fit.waverange));
legend_entry = cell(numel(fit.waverange), 1);
for i=1:numel(fit.waverange)
    [~, ind_w] = min(abs(handles.TA.WVec-fit.waverange(i)));
    currColor = ax4.ColorOrder(ax4.ColorOrderIndex,:);
    plot(ax4, handles.TA.TVec, fit.KineticFit(:,i)*1e3, 'Color', currColor, 'LineWidth', 2);
    scatter(ax4, handles.TA.TVec, handles.TA.deltaA(:,ind_w)*1e3, [], currColor, 'filled', 'SizeData', 10, 'MarkerFaceAlpha', 0.4);
    h(i) = plot(ax5, handles.TA.TVec, fit.KineticFit(:,i)*1e3, 'Color', currColor, 'LineWidth', 2);
    legend_entry{i} = strcat(string(round(handles.TA.WVec(ind_w))), ' nm');
    scatter(ax5, handles.TA.TVec, handles.TA.deltaA(:,ind_w)*1e3, [], currColor, 'filled', 'SizeData', 10, 'MarkerFaceAlpha', 0.4);
end
ylim(ax4, 'auto');
ylimmax = ylim(ax4);
ylimmax = ylimmax(2);
ylim(ax5, 'auto');
ylimmin = ylim(ax5);
ylimmin = ylimmin(1);
ylim([ylimmin, ylimmax]);
legend(ax5, h(:), legend_entry, 'Location', 'best');
legend('boxoff');
if strcmp(fs_or_ns, 'fsTA')
    xlabel('Time (ps)')
else
    xlabel('Time (\mus)')
end
ylabel(ax4, '\DeltaA (mOD)')
set([ax4,ax5],...
    'Box', 'on',...
    'LineWidth', 1.2,...
    'XMinorTick', 'on',...
    'YMinorTick', 'on',...
    'TickLength', [0.007 0.007],...
    'FontName', 'Helvetica',...
    'FontSize', 10);

ax6 = copyobj(handles.axes_spectrum, mainfig);
delete(ax6.Children(end));
spectrum_labels = cell(handles.poiList.Num_nodes,1);
for i=1:numel(spectrum_labels)
    pnt = handles.poiList.get_node(i);
    tInd = round(pnt.t);
    if handles.fsTA_radio.Value == 1
        spectrum_labels{i} = strcat(num2str(round(handles.TA.TVec(tInd),1)), ' ps');
    else
        spectrum_labels{i} = strcat(num2str(round(handles.TA.TVec(tInd),3)), ' \mus');
    end
end
legend(ax6, spectrum_labels, 'Location', 'best', 'Box', 'off');
set(ax6,...
    'Box', 'on',...
    'LineWidth', 1.2,...
    'XMinorTick', 'on',...
    'YMinorTick', 'on',...
    'TickLength', [0.007 0.007],...
    'FontName', 'Helvetica',...
    'FontSize', 10,...
    'Position', [0.0513 0.5656 0.4344 0.3901],...
    'Units', 'normalized');
            
end

