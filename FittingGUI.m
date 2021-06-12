function varargout = FittingGUI(varargin)
% FITTINGGUI MATLAB code for FittingGUI.fig
%      FITTINGGUI, by itself, creates a new FITTINGGUI or raises the existing
%      singleton*.
%
%      H = FITTINGGUI returns the handle to a new FITTINGGUI or the handle to
%      the existing singleton*.
%
%      FITTINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FITTINGGUI.M with the given input arguments.
%
%      FITTINGGUI('Property','Value',...) creates a new FITTINGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FittingGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FittingGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FittingGUI

% Last Modified by GUIDE v2.5 06-Apr-2021 17:45:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FittingGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @FittingGUI_OutputFcn, ...
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


% --- Executes just before FittingGUI is made visible.
function FittingGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FittingGUI (see VARARGIN)

% Choose default command line output for FittingGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% get data from main app
h = findobj('Tag','fig_surfxplorer');
handles.MainFigData = guidata(h);
TA = handles.MainFigData.TA;
handles.TA = TA;
guidata(hObject, handles);

table_data = {  't0',0, 0, false;...
                'IRF', 0.3, 0, false;...
                'tau1', 10, 0, false;};
handles.table_params.Data = table_data;
guidata(hObject, handles);

update_fitGUI_plots(handles,0);


% UIWAIT makes FittingGUI wait for user response (see UIRESUME)
% uiwait(handles.fig_fit_gui);


% --- Outputs from this function are returned to the command line.
function varargout = FittingGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in model_selector.
function model_selector_Callback(hObject, eventdata, handles)
% hObject    handle to model_selector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns model_selector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from model_selector

handles.export_fit.Enable = 'off';

contents = cellstr(get(hObject,'String'));
switch contents{get(hObject,'Value')}
    case {'A->Gnd', 'A->B'}
        table_data = {  't0',0,0,false;...
                'IRF', 0.3,0,false;...
                'tau1', 5,0,false;};
    case {'A->B->Gnd', 'A->B->C'}
        table_data = {  't0',0,0,false;...
                'IRF', 0.3,0,false;...
                'tau1', 5,0,false;...
                'tau2', 50,0,false;};
    case {'A->B->C->Gnd', 'A->B->C->D'}
        table_data = {  't0',0,0,false;...
                'IRF', 0.3,0,false;...
                'tau1', 5,0,false;...
                'tau2', 50,0,false;...
                'tau3', 500,0,false;};
    case 'A->B->C->D->Gnd'
        table_data = {  't0',0,0,false;...
                'IRF', 0.3,0,false;...
                'tau1', 5,0,false;...
                'tau2', 50,0,false;...
                'tau3', 500,0,false;...
                'tau4', 5000,0,false;};
    case 'Biexciton decay (1D) -> Gnd'
        table_data = {  't0',0,0,false;...
                'IRF', 0.3,0,false;...
                'OD(0)*k_A',1,0,false;};
    case 'Biexciton decay (1D) -> B'
        table_data = {  't0',0,0,false;...
                'IRF', 0.3,0,false;...
                'OD(0)*k_A',1,0,false;};
    case 'Biexciton decay (1D) -> Gnd with A->B'
        table_data = {  't0',0,0,false;...
                'IRF', 0.3,0,false;...
                'OD(0)*k_A', 1,0,false;
                'tau', 100,0,0};
    case 'Biexciton decay (1D) -> Gnd with A->Gnd'
        table_data = {  't0',0,0,false;...
            'IRF', 0.3,0,false;...
            'OD(0)*k_A', 1,0,false;
            'tau', 100,0,false};
end
handles.table_params.Data = table_data;
guidata(handles.fig_fit_gui, handles);
update_fitGUI_plots(handles,0);

% --- Executes during object creation, after setting all properties.
function model_selector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_selector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function waverange_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to waverange_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of waverange_textbox as text
%        str2double(get(hObject,'String')) returns contents of waverange_textbox as a double

update_fitGUI_plots(handles,0);
handles.menu_file_exportOrigin.Enable = 'off';

% --- Executes during object creation, after setting all properties.
function waverange_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to waverange_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in table_params.
function table_params_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table_params (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

update_fitGUI_plots(handles,0);
handles.menu_file_exportOrigin.Enable = 'off';


% --- Executes on button press in fit_button.
function fit_button_Callback(hObject, eventdata, handles)
% hObject    handle to fit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.fit = update_fitGUI_plots(handles,1);
handles.menu_file_exportOrigin.Enable = 'on';
guidata(handles.fig_fit_gui, handles);

% --- Executes on button press in export_fit.
function export_fit_Callback(hObject, eventdata, handles)
% hObject    handle to export_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % % contents = cellstr(get(handles.model_selector,'String'));
% % % switch contents{get(handles.model_selector,'Value')}
% % %     case {'A->Gnd', 'A->B'}
% % %         fitpar_labels = {'t0 (ps) = ','IRF (ps) = ','tau1 (ps) = '};
% % %     case {'A->B->Gnd', 'A->B->C'}
% % %         fitpar_labels = {'t0 (ps) = ','IRF (ps) = ','tau1 (ps) = ','tau2 (ps) = '};
% % %     case {'A->B->C->Gnd', 'A->B->C->D'}
% % %         fitpar_labels = {'t0 (ps) = ','IRF (ps) = ','tau1 (ps) = ','tau2 (ps) = ','tau3 (ps) = '};
% % %     case 'A->B->C->D->Gnd'
% % %         fitpar_labels = {'t0 (ps) = ','IRF (ps) = ','tau1 (ps) = ','tau2 (ps) = ','tau3 (ps) = ','tau4 (ps) = '};
% % %     case {'Biexciton decay (1D) -> Gnd', 'Biexciton decay (1D) -> B'}
% % %         fitpar_labels = {'t0 (ps) = ','IRF (ps) = ','kA*dOD(0) (ps^-0.5) = '};
% % %     case {'Biexciton decay (1D) -> Gnd with A->B', 'Biexciton decay (1D) -> Gnd with A->Gnd'}
% % %         fitpar_labels = {'t0 (ps) = ','IRF (ps) = ','kA*dOD(0) (ps^-0.5) = ','tau (ps) = '};
% % % end
% % % saveresults(handles.fit, fitpar_labels);


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_exportOrigin_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_exportOrigin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(handles.model_selector,'String'));
switch contents{get(handles.model_selector,'Value')}
    case {'A->Gnd', 'A->B'}
        fitpar_labels = {'t0 (ps) = ','IRF (ps) = ','tau1 (ps) = '};
    case {'A->B->Gnd', 'A->B->C'}
        fitpar_labels = {'t0 (ps) = ','IRF (ps) = ','tau1 (ps) = ','tau2 (ps) = '};
    case {'A->B->C->Gnd', 'A->B->C->D'}
        fitpar_labels = {'t0 (ps) = ','IRF (ps) = ','tau1 (ps) = ','tau2 (ps) = ','tau3 (ps) = '};
    case 'A->B->C->D->Gnd'
        fitpar_labels = {'t0 (ps) = ','IRF (ps) = ','tau1 (ps) = ','tau2 (ps) = ','tau3 (ps) = ','tau4 (ps) = '};
    case {'Biexciton decay (1D) -> Gnd', 'Biexciton decay (1D) -> B'}
        fitpar_labels = {'t0 (ps) = ','IRF (ps) = ','kA*dOD(0) (ps^-0.5) = '};
    case {'Biexciton decay (1D) -> Gnd with A->B', 'Biexciton decay (1D) -> Gnd with A->Gnd'}
        fitpar_labels = {'t0 (ps) = ','IRF (ps) = ','kA*dOD(0) (ps^-0.5) = ','tau (ps) = '};
end
exportFitToOrigin(handles.fit, fitpar_labels);

% --------------------------------------------------------------------
function menu_file_export_matlab_figures_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_export_matlab_figures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
exportPlots(handles,'Fitting');


% --------------------------------------------------------------------
function menu_file_loadModelLibrary_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_loadModelLibrary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected cell(s) is changed in table_params.
function table_params_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to table_params (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
