function varargout = SurfaceXplorerMatlab(varargin)
% SURFACEXPLORERMATLAB MATLAB code for SurfaceXplorerMatlab.fig
%      SURFACEXPLORERMATLAB, by itself, creates a new SURFACEXPLORERMATLAB or raises the existing
%      singleton*.
%
%      H = SURFACEXPLORERMATLAB returns the handle to a new SURFACEXPLORERMATLAB or the handle to
%      the existing singleton*.
%
%      SURFACEXPLORERMATLAB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SURFACEXPLORERMATLAB.M with the given input arguments.
%
%      SURFACEXPLORERMATLAB('Property','Value',...) creates a new SURFACEXPLORERMATLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SurfaceXplorerMatlab_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SurfaceXplorerMatlab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SurfaceXplorerMatlab

% Last Modified by GUIDE v2.5 05-Apr-2021 15:40:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SurfaceXplorerMatlab_OpeningFcn, ...
                   'gui_OutputFcn',  @SurfaceXplorerMatlab_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SurfaceXplorerMatlab is made visible.
function SurfaceXplorerMatlab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SurfaceXplorerMatlab (see VARARGIN)

% Choose default command line output for SurfaceXplorerMatlab
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% initialize dyanmic list of points-of-interest
poiList = POIList();
handles.poiList = poiList;

% initialize current working directory
handles.LastFolder = pwd;

% hides fitting text message
handles.fitting_text.Visible = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% just for debugging: uploading test.csv file for fast testing %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TA = TAsurface(strcat(handles.LastFolder, '/simulated_test.csv'));
handles.TA = TA;
updateGui(handles);
plotSurfAxes(handles);
handles = guidata(handles.fig_surfxplorer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

guidata(handles.fig_surfxplorer, handles);

% UIWAIT makes SurfaceXplorerMatlab wait for user response (see UIRESUME)
% uiwait(handles.fig_surfxplorer);


% --- Outputs from this function are returned to the command line.
function varargout = SurfaceXplorerMatlab_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function menu_file_open_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% open data file and read data
[FileName, PathName] = uigetfile({'*.csv';'*.txt'}, 'Select file to import', handles.LastFolder, 'MultiSelect', 'off');
if (PathName==0)
    return
end

% create TAsurface instance
TA = TAsurface(strcat(PathName, FileName));
handles.TA = TA;

% initialize points-of-interest list instance
poiList = POIList();
handles.poiList = poiList;

% initialize fit
handles.fit = cell(0);

% initialize GUI and plot surface
updateGui(handles);
plotSurfAxes(handles);

% save last folder for ease of importing new files from it
handles.LastFolder = PathName;

% save GUI struct
guidata(handles.fig_surfxplorer, handles);

% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_surface_Callback(hObject, eventdata, handles)
% hObject    handle to menu_surface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_surface_chirpCorrect_Callback(hObject, eventdata, handles)
% hObject    handle to menu_surface_chirpCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% open chirp correction GUI and wait for it to finish
chirpCorrect;
h = findobj('Tag','fig_chirpCorrect');
uiwait(h)

% chirp correct the data
t0Ind = getappdata(0,'t0Ind');
handles.TA.chirpCorrect(t0Ind);

% clear points-of-interest
button_clearPoints_Callback([],[],handles);

guidata(handles.fig_surfxplorer, handles);


% --- Executes on key press with focus on fig_surfxplorer and none of its controls.
function fig_surfxplorer_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to fig_surfxplorer (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% % % switch eventdata.Key
% % %     case 'g'
% % % end

%--------------------------------------------------------------------
function fig_surfxplorer_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to menu_addPointOfInterest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'TA')
    return;
end

% get cursor coordinates
p = get(handles.axes_surf, 'CurrentPoint');
wCursor = p(1,1);
tIndCursor = p(1,2);
% refresh 1d plots if cursor above 2d axes
if wCursor > min(handles.axes_surf.XLim) && ...
        wCursor < max(handles.axes_surf.XLim) && ...
        tIndCursor > min(handles.axes_surf.YLim) && ...
        tIndCursor < max(handles.axes_surf.YLim)
    refresh1DPlots(handles, wCursor, tIndCursor);
end


% --------------------------------------------------------------------
function menu_surface_subtractBackground_Callback(hObject, eventdata, handles)
% hObject    handle to menu_surface_subtractBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% open subtract background GUI and wait for it to finish
subtractBackground;
h = findobj('Tag','fig_subtractBackground');
uiwait(h)

% get background spectrum
background = getappdata(0,'background');
handles.TA.subtractBackground(background);

% clear points-of-interest
button_clearPoints_Callback([],[],handles);

guidata(handles.fig_surfxplorer, handles);

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function fig_surfxplorer_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to fig_surfxplorer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'TA')
    return;
end

% get cursor coordinates
p = get(handles.axes_surf, 'CurrentPoint');
wCursor = p(1,1);
tIndCursor = p(1,2);

% save 1d plots if cursor above 2d axes
if wCursor > min(handles.axes_surf.XLim) && ...
        wCursor < max(handles.axes_surf.XLim) && ...
        tIndCursor > min(handles.axes_surf.YLim) && ...
        tIndCursor < max(handles.axes_surf.YLim)
    
    handles.poiList.addPOI(wCursor, tIndCursor);
    plotSurfAxes(handles);
    % enable clear points button
    handles.button_clearPoints.Enable = 'on';
end


% --- Executes on button press in button_clearPoints.
function button_clearPoints_Callback(hObject, eventdata, handles)
% hObject    handle to button_clearPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% enable clear points button

handles.button_clearPoints.Enable = 'off';
handles.poiList.clearList();
plotSurfAxes(handles);
refresh1DPlots(handles, 1, 1);
% disable clear points button
handles.button_clearPoints.Enable = 'off';


% --------------------------------------------------------------------
function menu_fitting_svd_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting_svd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% open subtract background GUI and wait for it to finish

SVD;
h = findobj('Tag','fig_SVD');
uiwait(h)

% get principal component data
svdData = getappdata(0,'svdData');
handles.TA.SVDdata = svdData;

guidata(handles.fig_surfxplorer, handles);


% --------------------------------------------------------------------
function menu_surface_crop_Callback(hObject, eventdata, handles)
% hObject    handle to menu_surface_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_export_matlab_figures_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_export_matlab_figures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
exportPlots(handles,'SurfaceXplorerMatlab');

% --- Executes on button press in checkbox_kinetic_semilogx.
function checkbox_kinetic_semilogx_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_kinetic_semilogx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_kinetic_semilogx

% get cursor coordinates
p = get(handles.axes_surf, 'CurrentPoint');
wCursor = p(1,1);
tIndCursor = p(1,2);

refresh1DPlots(handles, wCursor, tIndCursor)


% --------------------------------------------------------------------
function menu_surface_quickCreatePOIs_Callback(hObject, eventdata, handles)
% hObject    handle to menu_surface_quickCreatePOIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% slice wavelength axis into equally-spaced wavelengths
w_inds = round(linspace(1, numel(handles.TA.WVec), 14));
poi_wavelengths = handles.TA.WVec(w_inds);

% create time indices for plotting and then add one for t<0 and for last
% time index
if handles.fsTA_radio.Value == 1
    poi_times = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000];
else
    poi_times = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10, 100, 200];
end
[~, t_ind] = min(abs(poi_times-handles.TA.TVec));
t_ind = [t_ind, numel(handles.TA.TVec)-2];
t_ind = unique(t_ind);
t_ind = [1, t_ind];

% remove existing points of interests
% clear all the already present points-of-interest
handles.button_clearPoints.Enable = 'off';
handles.poiList.clearList();

% add new points of interest
for i=1:numel(t_ind)
    handles.poiList.addPOI(poi_wavelengths(i), t_ind(i));
end

plotSurfAxes(handles);
refresh1DPlots(handles, 1, 1);

% re-enable clear points button
handles.button_clearPoints.Enable = 'on';

% --------------------------------------------------------------------
function menu_file_save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% prepare data struct for saving
deltaA = handles.TA.deltaA;
WVec = [0; handles.TA.WVec];
TVec = handles.TA.TVec;
TAsurface2save = [TVec, deltaA];
TAsurface2save = [WVec'; TAsurface2save]';

% save processed data in a new file
filter = {'*.csv'};
[FileName, PathName] = uiputfile('*.csv', 'Save processed data', handles.LastFolder);
if FileName==0
    return
end
csvwrite(strcat(PathName, FileName), TAsurface2save);

% --------------------------------------------------------------------
function menu_fitting_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_fitting_svd1stOrderFitting_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting_svd1stOrderFitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.TA.fit_to_model(3)


% --------------------------------------------------------------------
function menu_surface_t0_offset_Callback(hObject, eventdata, handles)
% hObject    handle to menu_surface_t0_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt = {'Enter t0 offset (ps for fsTA / ns for nsTA):'};
title = 'Time-zero offset';
dims = [1 50];
definput = {'0'};
answer = inputdlg(prompt,title,dims,definput);
if isempty(answer)
    return
end
if handles.fsTA_radio.Value == 1
    offset= str2num(answer{1});
else
    offset= str2num(answer{1})/1000;
end
handles.TA.t0_offset(offset);

plotSurfAxes(handles);


% --------------------------------------------------------------------
function menu_file_export_origin_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_export_origin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

deltaA = handles.TA.deltaA;
WVec = handles.TA.WVec;
TVec = handles.TA.TVec;

% get wavelength and times for exporting to origin
prompt = {'Enter wavelengths for kinetic traces (nm):','Enter times for spectral evolution (ps):'};
title = 'Input';
dims = [1 35];
definput = {'500 600 700 800','1 2 5 10 20 50 100 200 500 1000 2000 5000 7500'};
answer = inputdlg(prompt,title,dims,definput);
if isempty(answer)
    return
end
desired_wavelengths = str2num(answer{1});
desired_times = str2num(answer{2});

% find closest matching results for user selection
[~, wavelength_indices] = min(abs(WVec-desired_wavelengths));
selected_wavelengths = WVec(wavelength_indices)';

[~, time_indices] = min(abs(TVec-desired_times));
time_indices = [1, time_indices];
selected_times = TVec(time_indices);

% prepare data for saving
kinetics = [TVec, deltaA(:, wavelength_indices)];
spectra = [WVec, deltaA(time_indices, :)'];

% prepare headers: 3 header rows for: "Long Name", "Units", "Comments"
longName_kinetics = cell2mat(['Time', repmat({',\g(D)A'}, 1, numel(wavelength_indices))]);
units_kinetics = cell2mat(['ps', repmat({',OD'}, 1, numel(wavelength_indices))]);
comments_kinetics = ['0,', strjoin(strcat(strsplit(num2str(round(selected_wavelengths))), ' nm'),',')];

longName_spectra = cell2mat(['Wavelength', repmat({',\g(D)A'}, 1, numel(time_indices))]);
units_spectra = cell2mat(['nm', repmat({',OD'}, 1, numel(time_indices))]);
comments_spectra = ['0,', strjoin(strcat(strsplit(num2str(round(selected_times'))), ' ps'),',')];

% save as csv files with headers
[FileName, PathName] = uiputfile('*.csv', 'Choose master file name', handles.LastFolder);
if FileName == 0
    return
end
[PathName, FileName, ~] = fileparts(strcat(PathName, FileName));
% save kinetics
fid = fopen(strcat(PathName, '\', FileName, '_kinetics.csv'),'w'); 
fprintf(fid,'%s\n',longName_kinetics);
fprintf(fid,'%s\n',units_kinetics);
fprintf(fid,'%s\n',comments_kinetics);
fclose(fid);
dlmwrite(strcat(PathName, '\', FileName, '_kinetics.csv'), kinetics, '-append');
% save spectra
fid = fopen(strcat(PathName, '\', FileName, '_spectra.csv'),'w'); 
fprintf(fid,'%s\n',longName_spectra);
fprintf(fid,'%s\n',units_spectra);
fprintf(fid,'%s\n',comments_spectra);
fclose(fid);
dlmwrite(strcat(PathName, '\', FileName, '_spectra.csv'), spectra, '-append');

% --------------------------------------------------------------------
function menu_surface_pump_scatter_Callback(hObject, eventdata, handles)
% hObject    handle to menu_surface_pump_scatter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_fitting_sw1stOrderFitting_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting_sw1stOrderFitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SW_fit;
h = findobj('Tag','fig_SW_fit');
uiwait(h)

% % % handles.TA.SW_1st_order_fitting()


% --------------------------------------------------------------------
function menu_fitting_nlin1stOrderFitting_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting_nlin1stOrderFitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.TA.nlin_1st_order_fitting()


% --------------------------------------------------------------------
function menu_surface_quickCreatePOIs_spectrum_Callback(hObject, eventdata, handles)
% hObject    handle to menu_surface_quickCreatePOIs_spectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

time = handles.TA.TVec;
wavelength = handles.TA.WVec(round(end/2));

% find indices where time vector is greater than 1 ps
tg0_ind = time>1;

% create time indices for plotting and then add one for t<0 and for last
% time index
poi_times = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000];
[~, t_ind] = min(abs(poi_times-handles.TA.TVec));
t_ind = [t_ind, numel(handles.TA.TVec)-2];
t_ind = unique(t_ind);
t_ind = [1, t_ind];

% remove existing points of interests
% clear all the already present points-of-interest
handles.button_clearPoints.Enable = 'off';
handles.poiList.clearList();

% add new points of interest
for i=1:numel(t_ind)
    handles.poiList.addPOI(wavelength, t_ind(i));
end

plotSurfAxes(handles);

% re-enable clear points button
handles.button_clearPoints.Enable = 'on';

% --------------------------------------------------------------------
function menu_surface_quickCreatePOIs_kinetics_Callback(hObject, eventdata, handles)
% hObject    handle to menu_surface_quickCreatePOIs_kinetics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% find index in time vector closest to 1 ps
[~, t_ind] = min(abs(handles.TA.TVec-1));

% slice wavelength axis into equally-spaced wavelengths
w_inds = round(linspace(1, numel(handles.TA.WVec), 14));
wavelengths = handles.TA.WVec(w_inds);

% remove existing points of interests
% clear all the already present points-of-interest
handles.button_clearPoints.Enable = 'off';
handles.poiList.clearList();

% add new points of interest
for i = w_inds
    handles.poiList.addPOI(handles.TA.WVec(i), t_ind);
end

plotSurfAxes(handles);

% re-enable clear points button
handles.button_clearPoints.Enable = 'on';


% --------------------------------------------------------------------
function menu_surface_crop_keepView_Callback(hObject, eventdata, handles)
% hObject    handle to menu_surface_crop_keepView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.button_clearPoints.Enable = 'off';
handles.poiList.clearList();
wRange = handles.axes_surf.XLim;
tRange = round(handles.axes_surf.YLim);
if tRange(1)<1
    tRange(1)=1;
end
if tRange(2)>numel(handles.TA.TVec)
    tRange(2) = numel(handles.TA.TVec);
end
tRange = [handles.TA.TVec(tRange(1)), handles.TA.TVec(tRange(2))];
handles.TA.crop(wRange, tRange);
plotSurfAxes(handles);


% --------------------------------------------------------------------
function menu_surface_crop_deleteView_Callback(hObject, eventdata, handles)
% hObject    handle to menu_surface_crop_deleteView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.button_clearPoints.Enable = 'off';
handles.poiList.clearList();
wRange = handles.axes_surf.XLim;
handles.TA.deleteRegion(wRange);
plotSurfAxes(handles);

% --------------------------------------------------------------------
function menu_fitting_modelABC_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting_modelABC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_surface_quickCreatePOIs_Callback(handles.menu_surface_quickCreatePOIs, [], handles)

% get wavelengths to show in fit
prompt = {'Enter wavelengths to plot with fits (nm):'};
title = 'Input';
dims = [1 35];
definput = {''};
answer = inputdlg(prompt,title,dims,definput);
if isempty(answer)
    return
end
fit_wavelengths = str2num(answer{1});

handles.fitting_text.Visible = 'on';
drawnow();

data.spec = handles.TA.deltaA;
data.wavelength = handles.TA.WVec;
data.time = handles.TA.TVec;
dadt = @(a,k,t) [-k(1)*a(1);
                  k(1)*a(1)-k(2)*a(2);
                  k(2)*a(2)];
handles.model.kinmod = @(v) buildmodeldiff(data,dadt,v(1),v(2),[1/v(3) 1/v(4)],[1 0 0]);
handles.model.initialguess = [0 3 10 500];
handles.model.lb = [-2 0.1 0 0];
handles.model.ub = [2 5 2e4 2e4];
handles.fit = SWfit(data,fit_wavelengths,handles.model,'s');
if handles.fsTA_radio.Value == 1
    plotFit(handles, handles.fit, 'fsTA', 'A->B->C')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ps');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ps')});
else
    plotFit(handles, handles.fit, 'nsTA', 'A->B->C')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ns');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ns')});
end

handles.fitting_text.Visible = 'off';
guidata(handles.fig_surfxplorer, handles);


% --------------------------------------------------------------------
function menu_fitting_modelAB_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting_modelAB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_surface_quickCreatePOIs_Callback(handles.menu_surface_quickCreatePOIs, [], handles)

% get wavelengths to show in fit
prompt = {'Enter wavelengths to plot with fits (nm):'};
title = 'Input';
dims = [1 35];
definput = {''};
answer = inputdlg(prompt,title,dims,definput);
if isempty(answer)
    return
end
fit_wavelengths = str2num(answer{1});

handles.fitting_text.Visible = 'on';

data.spec = handles.TA.deltaA;
data.wavelength = handles.TA.WVec;
data.time = handles.TA.TVec;
dadt = @(a,k,t) [-k(1)*a(1);
                  k(1)*a(1);];
handles.model.kinmod = @(v) buildmodeldiff(data,dadt,v(1),v(2),[1/v(3)],[1 0]);
handles.model.initialguess = [0 3 50];
handles.model.lb = [-2 0.1 0];
handles.model.ub = [2 5 20000];
handles.fit = SWfit(data,fit_wavelengths,handles.model,'s');
if handles.fsTA_radio.Value == 1
    plotFit(handles, handles.fit, 'fsTA', 'A->B')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ps');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ps')});
else
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ns');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ns')});
end

handles.fitting_text.Visible = 'off';
guidata(handles.fig_surfxplorer, handles);


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_fitting_modelABG_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting_modelABG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_surface_quickCreatePOIs_Callback(handles.menu_surface_quickCreatePOIs, [], handles)

% get wavelengths to show in fit
prompt = {'Enter wavelengths to plot with fits (nm):'};
title = 'Input';
dims = [1 35];
definput = {''};
answer = inputdlg(prompt,title,dims,definput);
if isempty(answer)
    return
end
fit_wavelengths = str2num(answer{1});

handles.fitting_text.Visible = 'on';
drawnow();

data.spec = handles.TA.deltaA;
data.wavelength = handles.TA.WVec;
data.time = handles.TA.TVec;
dadt = @(a,k,t) [-k(1)*a(1);
                  k(1)*a(1)-k(2)*a(2);];
handles.model.kinmod = @(v) buildmodeldiff(data,dadt,v(1),v(2),[1/v(3) 1/v(4)],[1 0]);
handles.model.initialguess = [0 3 60 2000];
handles.model.lb = [-2 0.1 0 0];
handles.model.ub = [2 5 2e4 2e4];
handles.fit = SWfit(data,fit_wavelengths,handles.model,'s');
if handles.fsTA_radio.Value == 1
    plotFit(handles, handles.fit, 'fsTA', 'A->B->GND')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ps');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ps')});
else
    plotFit(handles, handles.fit, 'nsTA', 'A->B->GND')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ns');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ns')});
end

handles.fitting_text.Visible = 'off';
guidata(handles.fig_surfxplorer, handles);

% --------------------------------------------------------------------
function menu_fitting_modelABCG_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting_modelABCG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_surface_quickCreatePOIs_Callback(handles.menu_surface_quickCreatePOIs, [], handles)

% get wavelengths to show in fit
prompt = {'Enter wavelengths to plot with fits (nm):'};
title = 'Input';
dims = [1 35];
definput = {''};
answer = inputdlg(prompt,title,dims,definput);
if isempty(answer)
    return
end
fit_wavelengths = str2num(answer{1});

handles.fitting_text.Visible = 'on';
drawnow();

data.spec = handles.TA.deltaA;
data.wavelength = handles.TA.WVec;
data.time = handles.TA.TVec;
dadt = @(a,k,t) [-k(1)*a(1);
                  k(1)*a(1)-k(2)*a(2);
                  k(2)*a(2)-k(3)*a(3)];
handles.model.kinmod = @(v) buildmodeldiff(data,dadt,v(1),v(2),[1/v(3) 1/v(4) 1/v(5)],[1 0 0]);
handles.model.initialguess = [0 3 3 200 2000];
handles.model.lb = [-2 0.1 0 0 0];
handles.model.ub = [2 5 2e4 2e4 2e4];
handles.fit = SWfit(data,fit_wavelengths,handles.model,'s');
if handles.fsTA_radio.Value == 1
    plotFit(handles, handles.fit, 'fsTA', 'A->B->C->GND')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ps');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ps')});
else
    plotFit(handles, handles.fit, 'nsTA', 'A->B->C->GND')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ns');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ns')});
end

handles.fitting_text.Visible = 'off';
guidata(handles.fig_surfxplorer, handles);


% --------------------------------------------------------------------
function menu_fitting_model2AB_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting_model2AB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_surface_quickCreatePOIs_Callback(handles.menu_surface_quickCreatePOIs, [], handles)

% get wavelengths to show in fit
prompt = {'Enter wavelengths to plot with fits (nm):'};
title = 'Input';
dims = [1 35];
definput = {''};
answer = inputdlg(prompt,title,dims,definput);
if isempty(answer)
    return
end
fit_wavelengths = str2num(answer{1});

handles.fitting_text.Visible = 'on';
drawnow();

data.spec = handles.TA.deltaA;
data.wavelength = handles.TA.WVec;
data.time = handles.TA.TVec;
dadt = @(a,k,t) [-k(1)*a(1).^2;
                  k(1)*a(1).^2;];
handles.model.kinmod = @(v) buildmodeldiff(data,dadt,v(1),v(2),[1/v(3)],[1 0]);
handles.model.initialguess = [0 3 10];
handles.model.lb = [-2 0.1 0];
handles.model.ub = [2 5 2e4];
handles.fit = SWfit(data,fit_wavelengths,handles.model,'s');
if handles.fsTA_radio.Value == 1
    plotFit(handles, handles.fit, 'fsTA', '2A->B')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ps');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ps')});
else
    plotFit(handles, handles.fit, 'nsTA', '2A->B')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ns');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ns')});
end

handles.fitting_text.Visible = 'off';
guidata(handles.fig_surfxplorer, handles);


% --------------------------------------------------------------------
function menu_fitting_model2A2BC_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting_model2A2BC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_surface_quickCreatePOIs_Callback(handles.menu_surface_quickCreatePOIs, [], handles)

% get wavelengths to show in fit
prompt = {'Enter wavelengths to plot with fits (nm):'};
title = 'Input';
dims = [1 35];
definput = {''};
answer = inputdlg(prompt,title,dims,definput);
if isempty(answer)
    return
end
fit_wavelengths = str2num(answer{1});

handles.fitting_text.Visible = 'on';
drawnow();

data.spec = handles.TA.deltaA;
data.wavelength = handles.TA.WVec;
data.time = handles.TA.TVec;
dadt = @(a,k,t) [-k(1)*a(1);
                  k(1)*a(1)-k(2)*a(2)^2;
                  k(2)*a(2)^2;];
handles.model.kinmod = @(v) buildmodeldiff(data,dadt,v(1),v(2),[1/v(3) 1/v(4)],[1 0 0]);
handles.model.initialguess = [0 3 3 50];
handles.model.lb = [-2 0.1 0 0];
handles.model.ub = [2 5 2e4 2e4];
handles.fit = SWfit(data,fit_wavelengths,handles.model,'s');
if handles.fsTA_radio.Value == 1
    plotFit(handles, handles.fit, 'fsTA', 'A->B,2B->C')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ps');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ps')});
else
    plotFit(handles, handles.fit, 'nsTA', 'A->B,2B->C')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ns');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ns')});
end

handles.fitting_text.Visible = 'off';
guidata(handles.fig_surfxplorer, handles);


% --------------------------------------------------------------------
function menu_fitting_model2ABG_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting_model2ABG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_surface_quickCreatePOIs_Callback(handles.menu_surface_quickCreatePOIs, [], handles)

% get wavelengths to show in fit
prompt = {'Enter wavelengths to plot with fits (nm):'};
title = 'Input';
dims = [1 35];
definput = {''};
answer = inputdlg(prompt,title,dims,definput);
if isempty(answer)
    return
end
fit_wavelengths = str2num(answer{1});

handles.fitting_text.Visible = 'on';
drawnow();

data.spec = handles.TA.deltaA;
data.wavelength = handles.TA.WVec;
data.time = handles.TA.TVec;
dadt = @(a,k,t) [-k(1)*a(1).^2;
                  k(1)*a(1).^2-k(2)*a(2);];
handles.model.kinmod = @(v) buildmodeldiff(data,dadt,v(1),v(2),[1/v(3) 1/v(4)],[1 0]);
handles.model.initialguess = [0 3 3 2000];
handles.model.lb = [-2 0.1 0 0];
handles.model.ub = [2 5 2e4 2e4];
handles.fit = SWfit(data,fit_wavelengths,handles.model,'s');
if handles.fsTA_radio.Value == 1
    plotFit(handles, handles.fit, 'fsTA', '2A->B->GND')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ps');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ps')});
else
    plotFit(handles, handles.fit, 'nsTA', '2A->B->GND')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ns');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ns')});
end

handles.fitting_text.Visible = 'off';
guidata(handles.fig_surfxplorer, handles);


% --------------------------------------------------------------------
function menu_file_PlotInOrigin_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_PlotInOrigin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get wavelength and times for exporting to origin
prompt = {'Enter wavelengths for kinetic traces (nm):','Enter times for spectral evolution (ps):'};
title = 'Input';
dims = [1 35];
definput = {'500 600 700 800','1 2 5 10 20 50 100 200 500 1000 2000 5000 8000'};
answer = inputdlg(prompt,title,dims,definput);
if isempty(answer)
    return
end
desired_wavelengths = str2num(answer{1});
desired_times = str2num(answer{2});

% find closest matching results for user selection
[~, wavelength_indices] = min(abs(handles.TA.WVec-desired_wavelengths));
selected_wavelengths = handles.TA.WVec(wavelength_indices)';

[~, time_indices] = min(abs(handles.TA.TVec-desired_times));
time_indices = [1, time_indices];
selected_times = handles.TA.TVec(time_indices);

% initialize Origin server communication
originObj=actxserver('Origin.ApplicationSI');
invoke(originObj, 'Execute', 'doc -mc 1;');
invoke(originObj, 'IsModified', 'false');

% Load the custom project
invoke(originObj, 'Load', 'C:\Users\itais\Documents\MEGA\Matlab\SurfaceXplorerMatlab\CreatePlotInOrigin.opju');

% prepare data for plotting
spectra = [handles.TA.WVec, handles.TA.deltaA(time_indices, :)'*1000];
kinetics = [handles.TA.TVec, handles.TA.deltaA(:, wavelength_indices)*1000];

% Send this data over to the relevant worksheets
invoke(originObj, 'PutWorksheet', 'Spectra', spectra);
invoke(originObj, 'PutWorksheet', 'Kinetics', kinetics);

% set worksheet headers
spectra_worksheet = invoke(originObj, 'FindWorksheet', 'Spectra');
kinetics_worksheet = invoke(originObj, 'FindWorksheet', 'Kinetics');
spectra_worksheet_cols = invoke(spectra_worksheet, 'Columns');
kinetics_worksheet_cols = invoke(kinetics_worksheet, 'Columns');
for i=1:numel(selected_times)
    col = invoke(spectra_worksheet_cols, 'Item', uint8(i));
    if i==1
        invoke(col, 'Comments', strcat(num2str(round(selected_times(i))), ' ps'));
    else
        invoke(col, 'Comments', num2str(round(selected_times(i))));
    end
end
for i=1:numel(selected_wavelengths)
    col = invoke(kinetics_worksheet_cols, 'Item', uint8(i));
    invoke(col, 'Comments', strcat(num2str(round(selected_wavelengths(i))), ' nm'));
end

% rescale axes


% --------------------------------------------------------------------
function menu_surface_smooth_Callback(hObject, eventdata, handles)
% hObject    handle to menu_surface_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.TA.smooth_surface();
plotSurfAxes(handles);



function uplim_text_Callback(hObject, eventdata, handles)
% hObject    handle to uplim_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uplim_text as text
%        str2double(get(hObject,'String')) returns contents of uplim_text as a double

plotSurfAxes(handles);

% --- Executes during object creation, after setting all properties.
function uplim_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uplim_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lowlim_text_Callback(hObject, eventdata, handles)
% hObject    handle to lowlim_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowlim_text as text
%        str2double(get(hObject,'String')) returns contents of lowlim_text as a double

plotSurfAxes(handles);

% --- Executes during object creation, after setting all properties.
function lowlim_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowlim_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_surface_subtractRefSurf_Callback(hObject, eventdata, handles)
% hObject    handle to menu_surface_subtractRefSurf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% open data file and read data
[FileName, PathName] = uigetfile({'*.csv';'*.txt'}, 'Select file for reference surface', handles.LastFolder, 'MultiSelect', 'off');
if (PathName==0)
    return
end

% create TAsurface instance
ref = TAsurface(strcat(PathName, FileName));

% subtract the reference spectrum
handles.TA.subtractRef(ref);

plotSurfAxes(handles);


% --------------------------------------------------------------------
function menu_fitting_modelABCD_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting_modelABCD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_surface_quickCreatePOIs_Callback(handles.menu_surface_quickCreatePOIs, [], handles)

% get wavelengths to show in fit
prompt = {'Enter wavelengths to plot with fits (nm):'};
title = 'Input';
dims = [1 35];
definput = {''};
answer = inputdlg(prompt,title,dims,definput);
if isempty(answer)
    return
end
fit_wavelengths = str2num(answer{1});

handles.fitting_text.Visible = 'on';
drawnow();

data.spec = handles.TA.deltaA;
data.wavelength = handles.TA.WVec;
data.time = handles.TA.TVec;
dadt = @(a,k,t) [-k(1)*a(1);
                  k(1)*a(1)-k(2)*a(2);
                  k(2)*a(2)-k(3)*a(3);
                  k(3)*a(3)];
handles.model.kinmod = @(v) buildmodeldiff(data,dadt,v(1),v(2),[1/v(3) 1/v(4) 1/v(5)],[1 0 0 0]);
handles.model.initialguess = [0 3 10 500 1000];
handles.model.lb = [-2 0.1 0 0 0];
handles.model.ub = [2 5 2e4 2e4 2e4];
handles.fit = SWfit(data,fit_wavelengths,handles.model,'s');
if handles.fsTA_radio.Value == 1
    plotFit(handles, handles.fit, 'fsTA', 'A->B->C->D')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ps');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ps')});
else
    plotFit(handles, handles.fit, 'nsTA', 'A->B->C->D')
    msgbox({strcat('t0 = ',num2str(handles.fit.FitParams(1)),' ns');strcat('IRF = ',num2str(handles.fit.FitParams(2)),' ns')});
end

handles.fitting_text.Visible = 'off';
guidata(handles.fig_surfxplorer, handles);


% --------------------------------------------------------------------
function menu_fitting_opengui_Callback(hObject, eventdata, handles)
% hObject    handle to menu_fitting_opengui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FittingGUI;
h = findobj('Tag','fig_fit_gui');


% --------------------------------------------------------------------
function menu_surface_SVD_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to menu_surface_SVD_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% open chirp correction GUI and wait for it to finish

SVD;
h = findobj('Tag','fig_SVD');
uiwait(h);

% get background spectrum
svdData = getappdata(0,'svdData');
if ~isstruct(svdData)
    return;
else
    handles.TA.svdFilter(svdData);
end

% clear points-of-interest
button_clearPoints_Callback([],[],handles);

guidata(handles.fig_surfxplorer, handles);
