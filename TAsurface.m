classdef TAsurface < handle
    
    properties
        % general parameters
        Path;
        Header;
        
        % raw data
        deltaA;  %  transient absorption (logarithmically sampled in time)
        WVec;    %  wavelength vector
        TVec;    %  time vector
        
        % analysis data
        TVec_linsamp;   % linearly sampled time vector
        deltaA_linsamp; % linearly-sampled-in-time deltaA matrix
        IRF_FWHM;       %  FWHM of the instrument response function (assumed Gaussian)
        Populations;    %  population vector of the different species
        Spectra;        %  species associated spectra of the different species
        FitResiduals;   %  residual matrix of fit
        

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% PUBLIC CLASS METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = public)
        
        function self = TAsurface(filepath) 
            self.Header = [];
            self.Path = filepath;
            self.IRF_FWHM = 3;
            self.Populations = [];
            self.Spectra = [];
            self.FitResiduals = [];
            
            % open raw data
            rawdata=importdata(filepath);
            wavelength = rawdata(2:end,1);
            time = rawdata(1,2:end);
            if time(end)==time(end-1) % make sure time vector is monotonically increasing
                [~, maxind] = max(time);
                time = time(1:maxind);
                rawdata = rawdata(2:end,2:maxind+1);
            else
                rawdata = rawdata(2:end,2:end);
            end
            ind=~isnan(wavelength);
            wavelength=wavelength(ind);
            rawdata=rawdata(ind,:);
            
            % Interpolate the NaN points that might be present in the data
            ind=isinf(rawdata);
            rawdata(ind)=nan;
            
            % remove the lines that are all NaN from the edges
            ind=isnan(rawdata);
            fi=find(sum(ind,2)~=length(time),1,'first');
            li=find(sum(ind,2)~=length(time),1,'last');
            rawdata=rawdata(fi:li,:);
            wavelength=wavelength(fi:li);
            
            % Zero the NaN points on the edges
            ind=isnan(rawdata(1,:));
            rawdata(1,ind)=0;
            ind=isnan(rawdata(end,:));
            rawdata(end,ind)=0;
            
% % %             % Interpolate the NaN values
% % %             dat=rawdata(:);
% % %             dat(isnan(dat)) = interp1(find(~isnan(dat)), dat(~isnan(dat)), find(isnan(dat)), 'spline');
% % %             rawdata=reshape(dat,length(wavelength),length(time));
            rawdata = inpaint_nans(rawdata);
            
            % sometimes wavelength axis is reversed because of how the
            % spectrometer was aligned
            if wavelength(2)<wavelength(1)
                wavelength = flipud(wavelength);
                rawdata = flipud(rawdata);
            end
            
            % save as class properties
            self.deltaA = rawdata';
            self.WVec = wavelength;
            self.TVec = time';
            
% % %             % calculate linearly samples time vector
% % %             TVec_linsamp = min(self.TVec):self.TVec(2)-self.TVec(1):max(self.TVec);
% % %             self.TVec_linsamp = TVec_linsamp;
% % %             
% % %             % calculate linearly samples deltaA matrix
% % %             deltaA_linsamp = zeros(numel(TVec_linsamp), numel(self.deltaA(1,:)));
% % %             for i=1:numel(self.deltaA(1,:))
% % %                 deltaA_linsamp(:,i) = interp1(self.TVec, self.deltaA(:,i), TVec_linsamp);
% % %             end
% % %             self.deltaA_linsamp = deltaA_linsamp;
        end
        
        function subtractBackground(self, background)
            if numel(background)==numel(self.deltaA(1,:))
                self.deltaA = self.deltaA - repmat(background, numel(self.TVec), 1);
            else
                return
            end
        end
        
        function subtractRef(self, ref)
          
            [T,W] = meshgrid(ref.TVec, ref.WVec);
            [Tq,Wq] = meshgrid(self.TVec, self.WVec);
          
            % resample ref to fit size of current surface
            ref_deltaA_resamp = interp2(T, W, ref.deltaA', Tq, Wq)';
            
            self.deltaA = self.deltaA - ref_deltaA_resamp;
        end

        function chirpCorrect(self, shift_vector)
            if shift_vector(1) < 0
                return
            end
            TVec_linsamp_ = min(self.TVec):self.TVec(2)-self.TVec(1):max(self.TVec);
            deltaA_linsamp_ = zeros(numel(TVec_linsamp_), numel(self.deltaA(1,:)));
            for i=1:numel(self.deltaA(1,:))
                deltaA_linsamp_(:,i) = interp1(self.TVec, self.deltaA(:,i), TVec_linsamp_);
            end
            self.TVec_linsamp = TVec_linsamp_;
            self.deltaA_linsamp = deltaA_linsamp_;
            shift_vector = shift_vector-min(shift_vector);
            delta_size = round(max(shift_vector))+1; % the change in size of deltaA matrix after chirp correction
            wb = waitbar(0, 'Correcting chirp...');
            for i=1:numel(shift_vector)
                self.deltaA_linsamp(:,i) = fshift(self.deltaA_linsamp(:,i), -shift_vector(i));
                waitbar(i/(numel(shift_vector)+numel(self.deltaA_linsamp(1,:))), wb, 'Correcting chirp...');
            end
            for i=1:numel(self.deltaA_linsamp(1,:))
                self.deltaA(:,i) = interp1(self.TVec_linsamp, self.deltaA_linsamp(:,i), self.TVec);       
                waitbar((i+numel(shift_vector))/(numel(shift_vector)+numel(self.deltaA_linsamp(1,:))), wb, 'Saving corrected data...');
            end
            self.deltaA(isnan(self.deltaA)) = 0;
            close(wb)
        end
        
        function crop(self, wRange, tRange)
            wRange_ind = self.WVec >= wRange(1) & self.WVec <= wRange(2);
            tRange_ind = self.TVec >= tRange(1) & self.TVec <= tRange(2);
            self.deltaA = self.deltaA(tRange_ind, wRange_ind);
            self.WVec = self.WVec(wRange_ind);
            self.TVec = self.TVec(tRange_ind);
        end
        
        function deleteRegion(self, wRange)
            wRange_ind = self.WVec >= wRange(1) & self.WVec < wRange(2);
            self.deltaA(:, wRange_ind) = 0;
        end
        
        function t0_offset(self, offset)
            self.TVec = self.TVec-offset;
        end
        
        function smooth_surface(self)
            % smooth data with guassian filter
            self.deltaA = imgaussfilt(self.deltaA, 1);
        end
        
        function fit_to_model(self, model_name, fs_or_ns, fit_wavelengths, handles)
            
            % set initial (guess) value, number of model components
            % for SVD, and parameter vector to vary during minimization
            switch model_name
                case 'A->B'
                    k1 = 1/1;
                    params = [self.IRF_FWHM, k1];
                case 'A->B->C'
                    k1 = 1/1;
                    k2 = 1/10;
                    params = [self.IRF_FWHM, k1, k2];
                case 'A->B->GND'
                    k1 = 1/1;
                    k2 = 1/10;
                    params = [self.IRF_FWHM, k1, k2];
                case 'A->B->C->GND'
                    k1 = 1/1;
                    k2 = 1/10;
                    k3 = 1/100;
                    params = [self.IRF_FWHM, k1, k2, k3];
            end
            
            % set options for minimizer routine
            options=optimset('maxfunevals',3000, 'MaxIter', 3000,'TolFun',1e-10,'TolX',1e-10);
            
            % set residual calculating function
            resCalcFunc = @(x) self.calc_residual(x, model_name);
            
            % run minimization routine to fit the data
            [bestfit, ~, flag] = fminsearch(resCalcFunc, params, options);
            if flag<1
                disp('Bad initial guess or poor model');
                return
            end
            
            % calculate best fit parameters
            [~, self.Populations, self.Spectra, self.FitResiduals] = self.calc_residual(bestfit, model_name);
            self.Populations = self.Populations';
            
            % calculate fit errors
            [var, stdev, covMatrix] = self.fitStats(bestfit, model_name);
            
            % reconstructed deltaA from fit
            deltaArecon = self.Spectra*self.Populations';
                       
            % plot fit result
            mainfig = figure('Position', [192 72 980 591]);
            
            if strcmp(fs_or_ns, 'fsTA')
                ind = self.TVec < 1;
            else
                ind = self.TVec < 0.01;
            end
            ax1 = axes(mainfig, 'Position', [0.5552 0.0788 0.0695 0.3901], 'Units', 'normalized');
            plot(self.TVec, self.Populations/max(max(self.Populations)), 'LineWidth', 2);
            ylabel('Normalized model population')
            xlim([min(self.TVec) max(self.TVec(ind))])
            ax2 = axes(mainfig, 'Position', [0.6247 0.0788 0.3475 0.3901], 'Units', 'normalized');
            plot(self.TVec, self.Populations/max(max(self.Populations)), 'LineWidth', 2)
            xlim([min(self.TVec(~ind)), max(self.TVec(~ind))])
            set(gca,'Xscale','log','YTickLabel','')
            linkaxes([ax1 ax2],'y')
            if strcmp(fs_or_ns, 'fsTA')
                xlabel('Time (ps)')
            else
                xlabel('Time (\mus)')
            end
            legend('A','B','C', 'Location', 'best');
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
            plot(ax3, self.WVec, zeros(numel(self.WVec), 1), 'Color', [0.2,0.2,0.2]);
            ax3.ColorOrderIndex = 1;
            hold on
            h = plot(self.WVec, self.Spectra*1000, 'LineWidth', 2);
            xlabel('Wavelength (nm)')   
            ylabel('\DeltaA (mOD)')
            xlim([min(self.WVec), max(self.WVec)]) 
            switch model_name
                case 'A->B'
                    tauAB_txt = strcat('    \tau_{A->B}=', num2str(1./bestfit(end), '%.2f'), '\pm', num2str(stdev(end)/bestfit(end)^2, '%.2f'), time_unit);
                    legend(h, strcat('A ',tauAB_txt), 'B', 'Location', 'best');
                case 'A->B->C'
                    tauAB_txt = strcat('    \tau_{A->B}=', num2str(1./bestfit(end-1), '%.2f'), '\pm', num2str(stdev(end-1)/bestfit(end-1)^2, '%.2f'), time_unit);
                    tauBC_txt = strcat('    \tau_{B->C}=', num2str(1./bestfit(end), '%.2f'), '\pm', num2str(stdev(end)/bestfit(end)^2, '%.2f'), time_unit);
                    legend(h, strcat('A ',tauAB_txt), strcat('B ',tauBC_txt), 'C', 'Location', 'best');
                case 'A->B->GND'
                    tauAB_txt = strcat('    \tau_{A->B}=', num2str(1./bestfit(end-1), '%.2f'), '\pm', num2str(stdev(end-1)/bestfit(end-1)^2, '%.2f'), time_unit);
                    tauBG_txt = strcat('    \tau_{B->GND}=', num2str(1./bestfit(end), '%.2f'), '\pm', num2str(stdev(end)/bestfit(end)^2, '%.2f'), time_unit);
                    legend(h, strcat('A ',tauAB_txt), strcat('B ',tauBG_txt), 'Location', 'best');
                case 'A->B->C->GND'
                    tauAB_txt = strcat('    \tau_{A->B}=', num2str(1./bestfit(end-2), '%.2f'), '\pm', num2str(stdev(end-2)/bestfit(end-2)^2, '%.2f'), time_unit);
                    tauBC_txt = strcat('    \tau_{B->C}=', num2str(1./bestfit(end-1), '%.2f'), '\pm', num2str(stdev(end-1)/bestfit(end-1)^2, '%.2f'), time_unit);
                    tauCG_txt = strcat('    \tau_{C->GND}=', num2str(1./bestfit(end), '%.2f'), '\pm', num2str(stdev(end)/bestfit(end)^2, '%.2f'), time_unit);
                    legend(h, strcat('A ',tauAB_txt), strcat('B ',tauBC_txt), strcat('C ',tauCG_txt), 'Location', 'best');
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
                ind = self.TVec < 1;
            else
                ind = self.TVec < 0.01;
            end
            ax4 = axes(mainfig, 'Position', [0.0513 0.0788 0.0695 0.3901], 'Units', 'normalized'); 
            plot(ax4, self.TVec, zeros(numel(self.TVec), 1), 'Color', [0.2,0.2,0.2]);
            ax4.ColorOrderIndex = 1;
            hold on
            xlim([min(self.TVec) max(self.TVec(ind))]);
            ax5 = axes(mainfig, 'Position', [0.1208 0.0788 0.3649 0.3901], 'Units', 'normalized'); 
            plot(ax5, self.TVec, zeros(numel(self.TVec), 1), 'Color', [0.2,0.2,0.2]);
            ax5.ColorOrderIndex = 1;
            hold on
            xlim([min(self.TVec(~ind)), max(self.TVec(~ind))])
            set(gca,'Xscale','log','YTickLabel','')
            linkaxes([ax4 ax5],'y')
            h = zeros(1, numel(fit_wavelengths));
            legend_entry = cell(numel(fit_wavelengths), 1);
            for i=1:numel(fit_wavelengths)
                [~, ind_w] = min(abs(self.WVec-fit_wavelengths(i)));
                currColor = ax4.ColorOrder(ax4.ColorOrderIndex,:);
                plot(ax4, self.TVec, deltaArecon(ind_w,:)*1e3, 'Color', currColor, 'LineWidth', 2);
                scatter(ax4, self.TVec, self.deltaA(:,ind_w)'*1e3, [], currColor, 'filled', 'SizeData', 10, 'MarkerFaceAlpha', 0.4);
                h(i) = plot(ax5, self.TVec, deltaArecon(ind_w,:)*1e3, 'Color', currColor, 'LineWidth', 2);
                legend_entry{i} = strcat(string(round(self.WVec(ind_w))), ' nm');
                scatter(ax5, self.TVec, self.deltaA(:,ind_w)'*1e3, [], currColor, 'filled', 'SizeData', 10, 'MarkerFaceAlpha', 0.4);
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
            legend(ax6, spectrum_labels, 'NumColumns', 2, 'Location', 'best', 'Box', 'off');
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
                     
% % %             figure;
% % %             imagesc([min(self.WVec),max(self.WVec)], [1, numel(self.TVec)], self.FitResiduals*1e3);
% % %             ax = gca;
% % %             colormap(bluewhitered(128, ax.CLim));
% % %             title('Fit residuals (mOD)');
% % %             xlabel('Wavelength (nm)')
% % %             ylabel('Time (log scale)')
% % %             colorbar
% % %             ax.YTick = [];
% % %             ax.XLim = [min(self.WVec),max(self.WVec)];
% % %             ax.YDir = 'normal';
            
        end
        
        function svdFilter(self, svdData)
            self.deltaA = svdData.U*svdData.S*svdData.V';
        end
        
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% PRIVATE CLASS METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = public)
        
        function [residual, C, spec, residualMatrix] = calc_residual(self, params, model)
            
            % On the one hand we have the experimental A. On the other hand
            % we have A=spec*C (spectra times concentration). On each
            % iteration we vary the model parameters, and calculate the
            % resulting C, from it we calculate the trial TA spectra
            % (trialA=spec*C -> spec=AC'(CC')^(-1)) and then we try to
            % minimize (A-A_trial = A-spec*C)
            
            A = self.deltaA';
            
            % extend the time trace to remove end artifacts from gaussian
            % convolution
            dt = self.TVec(end)-self.TVec(end-1);
            time_extension = self.TVec(end)+(dt:dt:dt*numel(self.TVec));
            TVec_extended = [self.TVec; time_extension'];
            
            irfFWHM = params(1);
            switch model
                case 'A->B'
                    Kmatrix = [-params(2)      0; ...
                                params(2)      0];
                case 'A->B->C'
                    Kmatrix = [-params(2)      0              0; ...
                                params(2)     -params(3)      0; ...
                                0              params(3)      0];
                case 'A->B->GND'
                    Kmatrix = [-params(2)      0; ...
                                params(2)     -params(3)];
                case 'A->B->C->GND'
                    Kmatrix = [-params(2)      0              0; ...
                                params(2)     -params(3)      0; ...
                                0              params(3)     -params(4)];
            end
            
            % get the eigenvector matrix and the diagonal eigenvalue matrix
            % of the rate matrix
            [eigenVecMatrix, D] = eig(Kmatrix);
            
            % using knowledge on the populations at t=0 to find
            % constants for the initial value problem (IVP)
            switch model
                case 'A->B'
                    coeffs_of_IVP = eigenVecMatrix\[1;0];
                case 'A->B->C'
                    coeffs_of_IVP = eigenVecMatrix\[1;0;0];
                case 'A->B->GND'
                    coeffs_of_IVP = eigenVecMatrix\[1;0];
                case 'A->B->C->GND'
                    coeffs_of_IVP = eigenVecMatrix\[1;0;0];
            end
            
            % build solution for Populations at all times (without IRF
            % convolution)
            conc_vs_t = zeros(numel(Kmatrix(:,1)), numel(TVec_extended));
            for i=1:numel(TVec_extended)
                conc_vs_t(:,i) = eigenVecMatrix*expm(D*TVec_extended(i))*coeffs_of_IVP;
            end
            conc_vs_t(:, TVec_extended<0) = 0;
            
            % construct IRF
            b = 4*log(2)/irfFWHM^2;
            IRF = exp(-b.*self.TVec.^2);
            IRF = IRF/sum(IRF);
            
            % Perform the convolution and calculate C (see comment at
            % beginning of method code)
            for i=1:numel(Kmatrix(:,1))
                conv_output = conv(conc_vs_t(i,TVec_extended>0), IRF);
                C(i,:) = conv_output(1:numel(self.TVec));
            end
            
            % calculate the spectrum according to A=spec*C -> spec =
            % A*C'*(CC')^(-1) (remember C is not square)
            spec = A*C'/(C*C');
            
            % calculate the trial A
            trialA = spec*C;
            
            % calculate residual
            residualMatrix = (A-trialA)';
            
            % sum of the squares of residuals. For a non-square matrix M,
            % M'*M is square and its trace is the sum of squares of the
            % individual matrix components
            residual = trace(residualMatrix'*residualMatrix);
            
        end
        
        function [var, stdev, C] = fitStats(self, x, model)
            
            fit = self.Populations;
            var = (fit(:)'*fit(:))/(length(fit(:))-1);
                        
            % Calculate the jacobian matrix using a central difference approximation
            lx = length(x);
            J  = zeros(length(fit(:)),lx);
            for k = 1:lx
                if x(k)==0
                    dx = nthroot(eps,3);
                else
                    dx = x(k)*nthroot(eps,3);
                end
                xdp = x;
                xdm = x;
                xdp(k) = xdp(k)+dx;
                xdm(k) = xdm(k)-dx;
                [~, rd, ~, ~] = self.calc_residual(xdp, model);
                [~, rc, ~, ~] = self.calc_residual(xdm, model);
                J(:,k)=(rd(:)-rc(:))/(2*dx);
            end
            Hes = (J'*J);
            % Calculate the covariance matrix
            C = var*inv(Hes);
            % calculate the standard devation of the parameters
            stdev = sqrt(diag(C));
        end
        
    end
    
end
