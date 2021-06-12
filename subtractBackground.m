function varargout = subtractBackground(varargin)
% SUBTRACTBACKGROUND MATLAB code for subtractBackground.fig
%      SUBTRACTBACKGROUND, by itself, creates a new SUBTRACTBACKGROUND or raises the existing
%      singleton*.
%
%      H = SUBTRACTBACKGROUND returns the handle to a new SUBTRACTBACKGROUND or the handle to
%      the existing singleton*.
%
%      SUBTRACTBACKGROUND('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SUBTRACTBACKGROUND.M with the given input arguments.
%
%      SUBTRACTBACKGROUND('Property','Value',...) creates a new SUBTRACTBACKGROUND or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before subtractBackground_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to subtractBackground_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help subtractBackground

% Last Modified by GUIDE v2.5 03-Apr-2020 07:31:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @subtractBackground_OpeningFcn, ...
                   'gui_OutputFcn',  @subtractBackground_OutputFcn, ...
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


% --- Executes just before subtractBackground is made visible.
function subtractBackground_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to subtractBackground (see VARARGIN)

% Choose default command line output for subtractBackground
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

h = findobj('Tag','fig_surfxplorer');
handles.MainFigData = guidata(h);
TA = handles.MainFigData.TA;

plot(handles.axes_spectrum, TA.WVec, TA.deltaA(1, :)*1e3);
xlabel(handles.axes_spectrum, 'Wavelength (nm)')
ylabel(handles.axes_spectrum, '\DeltaA (mOD)')

handles.slider_numSpectra.Min = 1;
handles.slider_numSpectra.Max = numel(TA.TVec);
handles.slider_numSpectra.Value = 1;
handles.slider_numSpectra.SliderStep = 1/(numel(TA.TVec)-1)*ones(1,2);

% save mean spectrum
handles.meanSpectrum = TA.deltaA(1, :);

guidata(hObject, handles);

% UIWAIT makes subtractBackground wait for user response (see UIRESUME)
% uiwait(handles.fig_subtractBackground);


% --- Outputs from this function are returned to the command line.
function varargout = subtractBackground_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_numSpectra_Callback(hObject, eventdata, handles)
% hObject    handle to slider_numSpectra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

TA = handles.MainFigData.TA;

% set slider value to integer value
handles.slider_numSpectra.Value = round(handles.slider_numSpectra.Value);

% clear axes
cla(handles.axes_spectrum);

% plot number of plots indicated by slider
hold(handles.axes_spectrum, 'on')
for i=1:handles.slider_numSpectra.Value
    plot(handles.axes_spectrum, TA.WVec, TA.deltaA(i, :)*1e3);
end

% calculate mean spectrum and plot it
meanSpectrum = mean(TA.deltaA(1:handles.slider_numSpectra.Value, :), 1);
plot(handles.axes_spectrum, TA.WVec, meanSpectrum*1e3, 'LineWidth', 3, 'Color', 'k');
hold(handles.axes_spectrum, 'off')

% refresh number of spectra textbox
handles.text_numSpectra.String = num2str(handles.slider_numSpectra.Value);

% save mean spectrum
handles.meanSpectrum = meanSpectrum;
guidata(handles.fig_subtractBackground, handles);

% --- Executes during object creation, after setting all properties.
function slider_numSpectra_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_numSpectra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in button_accept.
function button_accept_Callback(hObject, eventdata, handles)
% hObject    handle to button_accept (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(0,'background',handles.meanSpectrum)
close(handles.fig_subtractBackground);


% --- Executes when user attempts to close fig_subtractBackground.
function fig_subtractBackground_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to fig_subtractBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
