function exportPlots(handles,guiType)

switch guiType
    
    case 'SurfaceXplorerMatlab'
        
        surface_figure = figure;
        kinetic_figure = figure('Position', [785 193 661 420]);
        spectrum_figure = figure('Position', [217 193 560 420]);

        surface_axes = copyobj(handles.axes_surf, surface_figure);
        kinetic_axes_lin = copyobj(handles.axes_kinetic_lin, kinetic_figure);
        kinetic_axes_log = copyobj(handles.axes_kinetic_log, kinetic_figure);
        spectrum_axes = copyobj(handles.axes_spectrum, spectrum_figure);

        % delete dynamic line (the one that is not stored in poiList but is updated
        % dynamically when the cursor hovers over the surface)
        delete(kinetic_axes_lin.Children(end));
        delete(kinetic_axes_log.Children(end));
        delete(spectrum_axes.Children(end));

        TA = handles.TA;
        wavelengths = TA.WVec;
        times = TA.TVec;

        % format legends
        kinetic_labels = cell(handles.poiList.Num_nodes,1);
        spectrum_labels = cell(handles.poiList.Num_nodes,1);
        for i=1:numel(kinetic_labels)
            pnt = handles.poiList.get_node(i);
            tInd = round(pnt.t);
            [~, wInd] = min(abs(TA.WVec - pnt.w));
            kinetic_labels{i} = strcat(num2str(round(wavelengths(wInd))), ' nm');
            if handles.fsTA_radio.Value == 1
                spectrum_labels{i} = strcat(num2str(round(times(tInd),1)), ' ps');
            else
                spectrum_labels{i} = strcat(num2str(round(times(tInd),3)), ' \mus');
            end
        end
        legend(kinetic_axes_log, kinetic_labels, 'NumColumns', 2, 'Location', 'best', 'Box', 'off');
        legend(spectrum_axes, spectrum_labels, 'NumColumns', 2, 'Location', 'best', 'Box', 'off');

        % format axes
        set(surface_axes,...
            'Box', 'on',...
            'LineWidth', 1.2,...
            'XMinorTick', 'on',...
            'YMinorTick', 'on',...
            'TickLength', [0.007 0.007],...
            'FontName', 'Lucida Bright',...
            'FontSize', 12,...
            'ActivePositionProperty', 'outerposition',...
            'Units', 'normalized',...
            'OuterPosition', [0 0 1 1],...
            'position',[0.1400 0.1200 0.7750 0.8150]);
        set(kinetic_axes_lin,...
            'Box', 'on',...
            'LineWidth', 1.2,...
            'XMinorTick', 'on',...
            'YMinorTick', 'on',...
            'TickLength', [0.007 0.007],...
            'FontName', 'Lucida Bright',...
            'FontSize', 12,...
            'ActivePositionProperty', 'outerposition',...
            'Units', 'normalized',...
            'OuterPosition', [0 0 1 1],...
            'position',[0.14 0.12857 0.17964 0.80643]);
        set(kinetic_axes_log,...
            'Box', 'on',...
            'LineWidth', 1.2,...
            'XMinorTick', 'on',...
            'YTickLabel', '',...
            'TickLength', [0.007 0.007],...
            'FontName', 'Lucida Bright',...
            'FontSize', 12,...
            'ActivePositionProperty', 'outerposition',...
            'Units', 'normalized',...
            'OuterPosition', [0 0 1 1],...
            'position',[0.32224 0.12857 0.59276 0.80643]);
        set(spectrum_axes,...
            'Box', 'on',...
            'LineWidth', 1.2,...
            'XMinorTick', 'on',...
            'YMinorTick', 'on',...
            'TickLength', [0.007 0.007],...
            'FontName', 'Lucida Bright',...
            'FontSize', 12,...
            'ActivePositionProperty', 'outerposition',...
            'Units', 'normalized',...
            'OuterPosition', [0 0 1 1],...
            'position',[0.1400 0.1200 0.7750 0.8150]);

        % % % cmap = bluewhitered(128, handles.axes_surf.CLim);
        % % % cmap = jet(100);
        % % % colormap(surface_figure, cmap)
        
    case 'Fitting'
        
        populations_fitting_figure = figure('Position', [274 170 560 420]);
        EAS_fitting_figure = figure('Position', [849 170 560 420]);
        kinetic_fitting_figure = figure('Position', [480 492 740 420]);
        
        kinetic_fitting_axes_lin = copyobj(handles.ax_lin_fit, kinetic_fitting_figure);
        kinetic_fitting_axes_log = copyobj(handles.ax_log_fit, kinetic_fitting_figure);
        populations_fitting_axes_lin = copyobj(handles.ax_populations_lin, populations_fitting_figure);
        populations_fitting_axes_log = copyobj(handles.ax_populations_log, populations_fitting_figure);
        EAS_fitting_axes = copyobj(handles.ax_EAS, EAS_fitting_figure);
        
        set(kinetic_fitting_axes_lin,...
            'Box', 'on',...
            'LineWidth', 1.2,...
            'XMinorTick', 'on',...
            'YMinorTick', 'on',...
            'TickLength', [0.007 0.007],...
            'FontName', 'Lucida Bright',...
            'FontSize', 12,...
            'ActivePositionProperty', 'outerposition',...
            'Units', 'normalized',...
            'OuterPosition', [0 0 1 1],...
            'position',[0.075922 0.1381 0.12798 0.7969]);
        set(kinetic_fitting_axes_log,...
            'Box', 'on',...
            'LineWidth', 1.2,...
            'XMinorTick', 'on',...
            'YMinorTick', 'on',...
            'TickLength', [0.007 0.007],...
            'FontName', 'Lucida Bright',...
            'FontSize', 12,...
            'ActivePositionProperty', 'outerposition',...
            'Units', 'normalized',...
            'OuterPosition', [0 0 1 1],...
            'position',[0.20716 0.1381 0.7538 0.7969]);
        set(populations_fitting_axes_lin,...
            'Box', 'on',...
            'LineWidth', 1.2,...
            'XMinorTick', 'on',...
            'YMinorTick', 'on',...
            'TickLength', [0.007 0.007],...
            'FontName', 'Lucida Bright',...
            'FontSize', 12,...
            'ActivePositionProperty', 'outerposition',...
            'Units', 'normalized',...
            'OuterPosition', [0 0 1 1],...
            'position',[0.15536 0.1381 0.15893 0.7969]);
        set(populations_fitting_axes_log,...
            'Box', 'on',...
            'LineWidth', 1.2,...
            'XMinorTick', 'on',...
            'YMinorTick', 'on',...
            'TickLength', [0.007 0.007],...
            'FontName', 'Lucida Bright',...
            'FontSize', 12,...
            'ActivePositionProperty', 'outerposition',...
            'Units', 'normalized',...
            'OuterPosition', [0 0 1 1],...
            'position',[0.31786 0.1381 0.6431 0.7969]);
        set(EAS_fitting_axes,...
            'Box', 'on',...
            'LineWidth', 1.2,...
            'XMinorTick', 'on',...
            'YMinorTick', 'on',...
            'TickLength', [0.007 0.007],...
            'FontName', 'Lucida Bright',...
            'FontSize', 12,...
            'ActivePositionProperty', 'outerposition',...
            'Units', 'normalized',...
            'OuterPosition', [0 0 1 1],...
            'position',[0.1125 0.13095 0.84846 0.80405]);
        
        % add legend
        waverange = str2num(handles.waverange_textbox.String);
        for i=1:2:numel(kinetic_fitting_axes_log.Children)
            p((i+1)/2) = kinetic_fitting_axes_log.Children(i);
        end
        legend(kinetic_fitting_axes_log, p, string(num2cell(waverange)));
end

end

