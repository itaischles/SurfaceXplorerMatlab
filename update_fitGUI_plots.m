function fit = update_fitGUI_plots(handles,fitflag)

TA = handles.TA;
userData = handles.table_params.Data;
experimentData.spec = TA.deltaA;
experimentData.wavelength = TA.WVec;
experimentData.time = TA.TVec;

waverange = str2num(handles.waverange_textbox.String);

% figure out which model to use
selectedMenuOption = cellstr(get(handles.model_selector,'String'));
modelName = selectedMenuOption{get(handles.model_selector,'Value')};
model = getModel(modelName,experimentData,userData);

% % % guidata(handles.fig_fit_gui, handles);

% prefit or full fit the data
if fitflag == 0 % pre-fit
    prefit = CheckModel(experimentData,model,waverange);
    fit = [];
else            % full fit
    fit = SWfit(experimentData,waverange,model,'s');
    if isfield(fit, 'FitParams') % fit succeeded
        handles.table_params.Data(:,2)=num2cell(fit.FitParams');
        handles.table_params.Data(:,3)=num2cell(fit.Std');
        fit.KineticFit = fit.KineticFit;
    else
        return;
    end
end

% clear and prepare axes
cla(handles.ax_lin_fit)
cla(handles.ax_log_fit)
cla(handles.ax_EAS)
cla(handles.ax_populations_lin)
cla(handles.ax_populations_log)
cla(handles.ax_residuals)
hold(handles.ax_lin_fit, 'on')
hold(handles.ax_log_fit, 'on')
hold(handles.ax_populations_lin, 'on')
hold(handles.ax_populations_log, 'on')
hold(handles.ax_EAS, 'on')

% prepare curve colors
cmap = colormap('lines');
colors = cmap(1:50,:);

% prepare lin/log break
linlog_break_time = 2; % (ps)
[~, break_ind] = min(abs(TA.TVec-linlog_break_time));

% plot kinetic curves with fits
w_ind = zeros(numel(waverange),1);
for i=1:numel(waverange)
    [~, w_ind(i)] = min(abs(TA.WVec-waverange(i)));
    scatter(handles.ax_lin_fit, TA.TVec(1:break_ind), TA.deltaA(1:break_ind, w_ind(i)), 15, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.25);
    scatter(handles.ax_log_fit, TA.TVec(break_ind+1:end), TA.deltaA(break_ind+1:end, w_ind(i)), 15, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.25);
    if fitflag == 0
        p(i) = plot(handles.ax_lin_fit, experimentData.time(1:break_ind), prefit(1:break_ind,i), 'Color', colors(i,:), 'LineWidth', 2);
        plot(handles.ax_log_fit, experimentData.time(break_ind+1:end), prefit(break_ind+1:end,i), 'Color', colors(i,:), 'LineWidth', 2);
    else
        p(i) = plot(handles.ax_lin_fit, experimentData.time(1:break_ind), fit.KineticFit(1:break_ind,i), 'Color', colors(i,:), 'LineWidth', 2);
        plot(handles.ax_log_fit, experimentData.time(break_ind+1:end), fit.KineticFit(break_ind+1:end,i), 'Color', colors(i,:), 'LineWidth', 2);
    end
    colorbar(handles.ax_residuals,'off')
end
if fitflag==1
    for i=1:numel(fit.Population(1,:))
        plot(handles.ax_populations_lin, experimentData.time(1:break_ind), fit.Population(1:break_ind,i)./max(fit.Population(:,i)), 'Color', colors(i,:), 'LineWidth', 2);
        plot(handles.ax_populations_log, experimentData.time(break_ind+1:end), fit.Population(break_ind+1:end,i)./max(fit.Population(:,i)), 'Color', colors(i,:), 'LineWidth', 2);
    end
    for i=1:numel(fit.SpectraFit(:,1))
        plot(handles.ax_EAS, experimentData.wavelength, fit.SpectraFit(i,:), 'Color', colors(i,:), 'LineWidth', 2);
    end
    ResMap = pcolor(handles.ax_residuals, TA.WVec, 1:1:numel(TA.TVec), fit.ResidualMap);
    set(ResMap, 'EdgeColor', 'none');
    shading(handles.ax_residuals,'interp');
    colormap(handles.ax_residuals, flipud(cbrewer('RdYlBu', 128, 'pchip')));
% % %     colormap(handles.ax_residuals, hot(128));
    colbar = colorbar(handles.ax_residuals);
    colbar.Label.String = '\DeltaA (OD)';
end

% format plots
handles.ax_lin_fit.YLim = [min(min(TA.deltaA(:, w_ind))),max(max(TA.deltaA(:, w_ind)))];
handles.ax_populations_lin.YLimMode = 'auto';
linkaxes([handles.ax_lin_fit, handles.ax_log_fit], 'y');
linkaxes([handles.ax_populations_lin, handles.ax_populations_log], 'y');
handles.ax_log_fit.XScale = 'log';
handles.ax_populations_log.XScale = 'log';
handles.ax_EAS.XLim = [min(experimentData.wavelength),max(experimentData.wavelength)];
handles.ax_log_fit.XLim = [linlog_break_time,max(TA.TVec)];
handles.ax_populations_log.XLim = [linlog_break_time,max(TA.TVec)];
handles.ax_log_fit.YTick = [];
handles.ax_populations_log.YTick = [];
handles.ax_lin_fit.FontSize = 10;
handles.ax_log_fit.FontSize = 10;
handles.ax_populations_lin.FontSize = 10;
handles.ax_populations_log.FontSize = 10;
handles.ax_EAS.FontSize = 10;
handles.ax_log_fit.XLabel.String = 'Time (ps)';
handles.ax_lin_fit.YLabel.String = '\DeltaA (OD)';
handles.ax_populations_log.XLabel.String = 'Time (ps)';
handles.ax_populations_lin.YLabel.String = 'Norm. populations';
handles.ax_EAS.XLabel.String = 'Wavelength (nm)';
handles.ax_EAS.YLabel.String = 'EAS';
handles.ax_residuals.XTick = '';
handles.ax_residuals.YTick = '';
handles.ax_residuals.XLabel.String = 'Wavelength';
handles.ax_residuals.YLabel.String = 'Time';
legend(handles.ax_log_fit, p, string(num2cell(waverange)));
if fitflag == 1
    legend(handles.ax_populations_log, model.species_names);
    legend(handles.ax_populations_log, 'boxoff');
    legend(handles.ax_EAS, model.species_names);
    legend(handles.ax_EAS, 'boxoff');
end

end

