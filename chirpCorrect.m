function varargout = chirpCorrect(varargin)
% CHIRPCORRECT MATLAB code for chirpCorrect.fig
%      CHIRPCORRECT, by itself, creates a new CHIRPCORRECT or raises the existing
%      singleton*.
%
%      H = CHIRPCORRECT returns the handle to a new CHIRPCORRECT or the handle to
%      the existing singleton*.
%
%      CHIRPCORRECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHIRPCORRECT.M with the given input arguments.
%
%      CHIRPCORRECT('Property','Value',...) creates a new CHIRPCORRECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before chirpCorrect_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to chirpCorrect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu_file.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help chirpCorrect

% Last Modified by GUIDE v2.5 02-Dec-2019 10:43:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @chirpCorrect_OpeningFcn, ...
                   'gui_OutputFcn',  @chirpCorrect_OutputFcn, ...
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


% --- Executes just before chirpCorrect is made visible.
function chirpCorrect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to chirpCorrect (see VARARGIN)

% Choose default command line output for chirpCorrect
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

handles.button_acceptFit.Enable = 'off';

handles.uitable_pts.Data = cell(0,3);

h = findobj('Tag','fig_surfxplorer');
handles.MainFigData = guidata(h);
TA = handles.MainFigData.TA;

% linear time sampling
% % % dtmin = min(diff(TA.TVec));
% % % TVec_linsamp = min(TA.TVec):dtmin:max(TA.TVec);
TVec_linsamp = min(TA.TVec):TA.TVec(2)-TA.TVec(1):max(TA.TVec);
deltaA_linsamp = zeros(numel(TVec_linsamp), numel(TA.deltaA(1,:)));

for i=1:numel(TA.deltaA(1,:))
    deltaA_linsamp(:,i) = interp1(TA.TVec, TA.deltaA(:,i), TVec_linsamp);
end

handles.TVec_linsamp = TVec_linsamp;
handles.deltaA_linsamp = deltaA_linsamp;
handles.MainFigData = guidata(h);

scaled_TA = atan(TA.deltaA*1000);
TAsurf = pcolor(handles.axes_surf, TA.WVec, 1:1:numel(TA.TVec), scaled_TA);
set(TAsurf, 'EdgeColor', 'none');
shading interp

% % % imagesc(handles.axes_surf, [min(TA.WVec),max(TA.WVec)], [1, numel(TVec_linsamp)], deltaA_linsamp);
handles.axes_surf.XLim = [min(TA.WVec),max(TA.WVec)];
[~, t_ind] = min(abs(TVec_linsamp - 10));
handles.axes_surf.YLim = [1,t_ind];
handles.axes_surf.YDir = 'normal';
xlabel(handles.axes_surf, 'Wavelength (nm)');
ylabel(handles.axes_surf, '');
% % % cmap = bluewhitered(128, handles.axes_surf.CLim);
% % % colormap(handles.fig_chirpCorrect, cmap)
cmap = flipud(cbrewer('RdBu', 128, 'pchip'));
colormap(handles.axes_surf, cmap)

guidata(handles.fig_chirpCorrect, handles);

% UIWAIT makes chirpCorrect wait for user response (see UIRESUME)
% uiwait(handles.fig_chirpCorrect);


% --- Outputs from this function are returned to the command line.
function varargout = chirpCorrect_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_addPt.
function button_addPt_Callback(hObject, eventdata, handles)
% hObject    handle to button_addPt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hold on

TA = handles.MainFigData.TA;

axes(handles.axes_surf);
[w_cursor,t0_cursor] = ginput(1);
handles.uitable_pts.Data{end+1,1} = w_cursor;
handles.uitable_pts.Data{end,2} = handles.TVec_linsamp(round(t0_cursor));
handles.uitable_pts.Data{end,3} = t0_cursor;

handles.plotPts = plot(handles.axes_surf, w_cursor, t0_cursor, 'Marker', 'o', 'Color', 'k', 'LineStyle', 'none');

% fit if enough pts are on plot
numRows = size(handles.uitable_pts.Data);
numRows = numRows(1);
if numRows >= 4
    handles.button_acceptFit.Enable = 'on';
    % delete previous fit if exists
    if isfield(handles, 'fitplot')
        delete(handles.fitplot)
    end
    % prepare data for fitting
    pts = cell2mat(handles.uitable_pts.Data);
    w = TA.WVec;
    x = pts(:,1);
    y = pts(:,3);
% % %     % prepare fitting model and parameters
% % %     ft = fittype( 'a*x.^b+c', 'independent', 'x', 'dependent', 'y' );
% % %     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% % %     opts.Display = 'Off';
% % %     opts.StartPoint = [0, 1, 0];
% % %     pts = cell2mat(handles.uitable_pts.Data);
% % %     % perform fit
% % %     fitresult = fit( x, y, ft, opts );
% % %     handles.fitresult = fitresult;
% % %     a = fitresult.a;
% % %     b = fitresult.b;
% % %     c = fitresult.c;
% % %     % plot resulting fit
% % %     t0Ind = a*(w/min(w)).^b+c;

    t0Ind = interp1(x, y, w, 'pchip');


    handles.fitplot = plot(handles.axes_surf, w, t0Ind, 'k');
    handles.t0Ind = t0Ind;
    guidata(handles.fig_chirpCorrect, handles)
end


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on fig_chirpCorrect and none of its controls.
function fig_chirpCorrect_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to fig_chirpCorrect (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

switch eventdata.Key
    case 'a'
        button_addPt_Callback([], [], handles)
end


% --- Executes on button press in button_clearPts.
function button_clearPts_Callback(hObject, eventdata, handles)
% hObject    handle to button_clearPts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hold off

TA = handles.MainFigData.TA;
TVec_linsamp = handles.TVec_linsamp;
deltaA_linsamp = handles.deltaA_linsamp;

% % % imagesc(handles.axes_surf, [min(TA.WVec),max(TA.WVec)], [1, numel(TVec_linsamp)], deltaA_linsamp);
scaled_TA = atan(TA.deltaA*1000);
TAsurf = pcolor(handles.axes_surf, TA.WVec, 1:1:numel(TA.TVec), scaled_TA);
set(TAsurf, 'EdgeColor', 'none');
shading interp

cmap = flipud(cbrewer('div', 'RdBu', 128, 'pchip'));
colormap(handles.axes_surf, cmap)
handles.axes_surf.YTick = [];
handles.axes_surf.XLim = [min(TA.WVec),max(TA.WVec)];
[~, t_ind] = min(abs(TVec_linsamp - 10));
handles.axes_surf.YLim = [1,t_ind];
handles.axes_surf.YDir = 'normal';
xlabel(handles.axes_surf, 'Wavelength (nm)');
ylabel(handles.axes_surf, 'Time (ps)');

handles.uitable_pts.Data = cell(0,3);

handles.button_acceptFit.Enable = 'off';


% --- Executes on button press in button_acceptFit.
function button_acceptFit_Callback(hObject, eventdata, handles)
% hObject    handle to button_acceptFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(0,'t0Ind',handles.t0Ind)
delete(handles.fig_chirpCorrect);


% --- Executes when user attempts to close fig_chirpCorrect.
function fig_chirpCorrect_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to fig_chirpCorrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
setappdata(0,'t0Ind',-1)
delete(hObject);
