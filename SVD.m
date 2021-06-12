function varargout = SVD(varargin)
% SVD MATLAB code for SVD.fig
%      SVD, by itself, creates a new SVD or raises the existing
%      singleton*.
%
%      H = SVD returns the handle to a new SVD or the handle to
%      the existing singleton*.
%
%      SVD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SVD.M with the given input arguments.
%
%      SVD('Property','Value',...) creates a new SVD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SVD_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SVD_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SVD

% Last Modified by GUIDE v2.5 05-Apr-2021 15:49:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SVD_OpeningFcn, ...
                   'gui_OutputFcn',  @SVD_OutputFcn, ...
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


% --- Executes just before SVD is made visible.
function SVD_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SVD (see VARARGIN)

% Choose default command line output for SVD
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

h = findobj('Tag','fig_surfxplorer');
handles.MainFigData = guidata(h);
TA = handles.MainFigData.TA;

% get desired number of principal components
N = str2num(handles.edit_numPrncplComp.String);

% do the SVD algorithm
[Ureduced, Sreduced, Vreduced] = svds(TA.deltaA, N);
Arecon = Ureduced*Sreduced*Vreduced';
Sdiag = diag(Sreduced,0);

% calculate delta(deltaA)^2
ddA = (Arecon - TA.deltaA).^2;

% plot delta(deltaA)^2
imagesc(handles.axes_ddA, ddA);
handles.axes_ddA.YTick = [];
handles.axes_ddA.XTick = [];
handles.axes_ddA.YDir = 'normal';
cmap = bluewhitered(128, handles.axes_ddA.CLim);
colormap(handles.fig_SVD, cmap)

% plot principal components
cla(handles.axes_U_lin)
cla(handles.axes_U_log)
cla(handles.axes_V)
clear handles.LegendU
clear handles.LegendV
hold(handles.axes_U_lin, 'on')
hold(handles.axes_U_log, 'on')
hold(handles.axes_V, 'on')
for i=1:N
    plot(handles.axes_U_lin, TA.TVec(TA.TVec<=2), Ureduced(1:numel(TA.TVec(TA.TVec<=2)),i)*Sdiag(i));
    plot(handles.axes_U_log, TA.TVec(TA.TVec>2), Ureduced(numel(TA.TVec(TA.TVec<=2))+1:end,i)*Sdiag(i));
    plot(handles.axes_V, TA.WVec, Vreduced(:,i)*Sdiag(i));
end
handles.LegendU = legend(handles.axes_U_log, num2str(Sdiag));
handles.LegendV = legend(handles.axes_V, num2str(Sdiag));
hold(handles.axes_U_lin, 'off')
hold(handles.axes_U_log, 'off')
hold(handles.axes_V, 'off')
handles.axes_U_log.YTick = [];
handles.axes_U_log.XScale = 'log';

% store values
handles.U = Ureduced;
handles.S = Sreduced;
handles.V = Vreduced;

guidata(hObject, handles);


% UIWAIT makes SVD wait for user response (see UIRESUME)
% uiwait(handles.fig_SVD);


% --- Outputs from this function are returned to the command line.
function varargout = SVD_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_numPrncplComp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numPrncplComp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TA = handles.MainFigData.TA;

% get desired number of principal components
N = round(str2num(handles.edit_numPrncplComp.String));

if N<1
    return;
end

% do the SVD algorithm
[Ureduced, Sreduced, Vreduced] = svds(TA.deltaA, N);
Arecon = Ureduced*Sreduced*Vreduced';
Sdiag = diag(Sreduced,0);

% calculate delta(deltaA)^2
ddA = (Arecon - TA.deltaA).^2;

% plot delta(deltaA)^2
imagesc(handles.axes_ddA, ddA);
handles.axes_ddA.YTick = [];
handles.axes_ddA.XTick = [];
handles.axes_ddA.YDir = 'normal';
cmap = bluewhitered(128, handles.axes_ddA.CLim);
colormap(handles.fig_SVD, cmap)

% plot principal components
cla(handles.axes_U_lin)
cla(handles.axes_U_log)
cla(handles.axes_V)
clear handles.LegendU
clear handles.LegendV
hold(handles.axes_U_lin, 'on')
hold(handles.axes_U_log, 'on')
hold(handles.axes_V, 'on')
for i=1:N
    plot(handles.axes_U_lin, TA.TVec(TA.TVec<=2), Ureduced(1:numel(TA.TVec(TA.TVec<=2)),i)*Sdiag(i));
    plot(handles.axes_U_log, TA.TVec(TA.TVec>2), Ureduced(numel(TA.TVec(TA.TVec<=2))+1:end,i)*Sdiag(i));
    plot(handles.axes_V, TA.WVec, Vreduced(:,i)*Sdiag(i));
end
handles.LegendU = legend(handles.axes_U_log, num2str(Sdiag));
handles.LegendV = legend(handles.axes_V, num2str(Sdiag));
hold(handles.axes_U_lin, 'off')
hold(handles.axes_U_log, 'off')
hold(handles.axes_V, 'off')
handles.axes_U_log.YTick = [];
handles.axes_U_log.XScale = 'log';

% store values
handles.U = Ureduced;
handles.S = Sreduced;
handles.V = Vreduced;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_numPrncplComp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numPrncplComp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_SVD_filter.
function button_SVD_filter_Callback(hObject, eventdata, handles)
% hObject    handle to button_SVD_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
svdData.U = handles.U; % principal kinetic components
svdData.S = handles.S; % singular values
svdData.V = handles.V; % principal spectral components
setappdata(0,'svdData',svdData)
close(handles.fig_SVD);


% --- Executes when user attempts to close fig_SVD.
function fig_SVD_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to fig_SVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.fig_SVD);


% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to button_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
svdData = 0;
setappdata(0,'svdData',svdData)
close(handles.fig_SVD);