function varargout = GUI_Transit_Time_3_signals(varargin)
% GUI_TRANSIT_TIME_3_SIGNALS MATLAB code for GUI_Transit_Time_3_signals.fig
%      GUI_TRANSIT_TIME_3_SIGNALS, by itself, creates a new GUI_TRANSIT_TIME_3_SIGNALS or raises the existing
%      singleton*.
%
%      H = GUI_TRANSIT_TIME_3_SIGNALS returns the handle to a new GUI_TRANSIT_TIME_3_SIGNALS or the handle to
%      the existing singleton*.
%
%      GUI_TRANSIT_TIME_3_SIGNALS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_TRANSIT_TIME_3_SIGNALS.M with the given input arguments.
%
%      GUI_TRANSIT_TIME_3_SIGNALS('Property','Value',...) creates a new GUI_TRANSIT_TIME_3_SIGNALS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Transit_Time_3_signals_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Transit_Time_3_signals_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Transit_Time_3_signals

% Last Modified by GUIDE v2.5 14-Jan-2021 17:24:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Transit_Time_3_signals_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Transit_Time_3_signals_OutputFcn, ...
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


% --- Executes just before GUI_Transit_Time_3_signals is made visible.
function GUI_Transit_Time_3_signals_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Transit_Time_3_signals (see VARARGIN)

global PP Amin_asc Amax_asc Amin_desc Amax_desc Amin_dia Amax_dia areamm2...
    Dis_Asc Dis_Desc Dis_Dia EDV EF  DBP SBP Pes ENd_est Ees SV tNd ENd_avg
PP = []; Amin_asc = []; Amax_asc = []; Amin_desc = []; Amax_desc = []; Amin_dia = [];...
    Amax_dia = []; areamm2 = []; Dis_Asc = []; Dis_Desc = []; Dis_Dia = []; EDV = [];...
    EF = []; DBP = []; SBP = []; Pes = []; ENd_est = []; Ees = []; SV = []; tNd = []; ENd_avg = [];

% Choose default command line output for GUI_Transit_Time_3_signals
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Transit_Time_3_signals wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_Transit_Time_3_signals_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
close 

% --- Executes on button press in Save_button.
function Save_button_Callback(hObject, eventdata, handles)
global transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia Dis_Asc Dis_Desc Dis_Dia EF1 SV HR CO info Q
info = getappdata(0,'info')
setappdata(0,'transit_time_asc_desc',transit_time_asc_desc)
setappdata(0,'transit_time_desc_dia',transit_time_desc_dia)
setappdata(0,'transit_time_asc_dia',transit_time_asc_dia)
%- Generate spreadsheet
Generate_spreadsheet_button_Callback

uiresume
% Update handles structure
guidata(hObject, handles);

% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
global transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia PP areamm2 Dis_Asc Dis_Desc Dis_Dia EF1 SV HR CO tNd
flow = getappdata(0,'flow');
freq = getappdata(0,'freq');
areamm2 = getappdata(0,'area');
waveform = 2; %- velocity/flow waves
continuous = 0; %- Single wave
show = 1; %- Display wave in handles.axes
if max(flow(:,2))<abs(min(flow(:,2))) %- Check descending flow wave is positive
    flow(:,2) = -flow(:,2);
end
if max(flow(:,3))<abs(min(flow(:,3))) %- Check diaphragm flow wave is positive
    flow(:,3) = -flow(:,3);
end
%- Make sure flows end at 0
for k=1:size(flow,2)
   temp = flow(:,k);
   d = temp(end);
   temp = temp-d;
   flow(:,k)=temp;
end
    
ymax_array = max(flow); ymax = max(ymax_array);
ymin_array = min(flow); ymin = min(ymin_array);
switch(get(eventdata.NewValue,'Tag'));
    case 'F2F_button'
        arrayfun(@cla,findall(0,'type','axes')),set(handles.str_TT_asc_desc, 'String', []),set(handles.str_TT_desc_dia, 'String', []),set(handles.str_TT_asc_dia, 'String', []);
        algorithm = 1;
        %- Ascending-Descending
        flow_asc_desc = [flow(:,1) flow(:,2)];
        flag_plot = 'asc_desc';
        transit_time_asc_desc = TTAlgorithm(flow_asc_desc,freq,algorithm,waveform,continuous,show,handles.axes_asc_desc,flag_plot);
        TT_asc_desc_str = sprintf('%.5f', transit_time_asc_desc);
        set(handles.str_TT_asc_desc, 'String', TT_asc_desc_str);
        axes(handles.axes_asc_desc), ylim([ymin*1.2 ymax*1.2]), 
        %- Descending-Diaphragm
        flow_desc_dia = [flow(:,2) flow(:,3)];
        flag_plot = 'desc_dia';
        transit_time_desc_dia = TTAlgorithm(flow_desc_dia,freq,algorithm,waveform,continuous,show,handles.axes_desc_dia,flag_plot);
        TT_desc_dia_str = sprintf('%.5f', transit_time_desc_dia);
        set(handles.str_TT_desc_dia, 'String', TT_desc_dia_str);
        axes(handles.axes_desc_dia), ylim([ymin*1.2 ymax*1.2])
        %- Ascending-Diaphragm
        flow_asc_dia = [flow(:,1) flow(:,3)];
        flag_plot = 'asc_dia';
        transit_time_asc_dia = TTAlgorithm(flow_asc_dia,freq,algorithm,waveform,continuous,show,handles.axes_asc_dia,flag_plot);
        transit_time_asc_dia = transit_time_asc_desc + transit_time_desc_dia;
        TT_asc_dia_str = sprintf('%.5f', transit_time_asc_dia);
        set(handles.str_TT_asc_dia, 'String', TT_asc_dia_str);
        axes(handles.axes_asc_dia), ylim([ymin*1.2 ymax*1.2])
    case 'F2F_radius_button'
        arrayfun(@cla,findall(0,'type','axes')),set(handles.str_TT_asc_desc, 'String', []),set(handles.str_TT_desc_dia, 'String', []),set(handles.str_TT_asc_dia, 'String', []);
        algorithm = 2;
        %- Ascending-Descending
        flow_asc_desc = [flow(:,1) flow(:,2)];
        flag_plot = 'asc_desc';
        transit_time_asc_desc = TTAlgorithm(flow_asc_desc,freq,algorithm,waveform,continuous,show,handles.axes_asc_desc,flag_plot);
        TT_asc_desc_str = sprintf('%.5f', transit_time_asc_desc);
        set(handles.str_TT_asc_desc, 'String', TT_asc_desc_str);
        axes(handles.axes_asc_desc), ylim([ymin*1.2 ymax*1.2]), 
        %- Descending-Diaphragm
        flow_desc_dia = [flow(:,2) flow(:,3)];
        flag_plot = 'desc_dia';
        transit_time_desc_dia = TTAlgorithm(flow_desc_dia,freq,algorithm,waveform,continuous,show,handles.axes_desc_dia,flag_plot);
        TT_desc_dia_str = sprintf('%.5f', transit_time_desc_dia);
        set(handles.str_TT_desc_dia, 'String', TT_desc_dia_str);
        axes(handles.axes_desc_dia), ylim([ymin*1.2 ymax*1.2])
        %- Ascending-Diaphragm
        flow_asc_dia = [flow(:,1) flow(:,3)];
        flag_plot = 'asc_dia';
        transit_time_asc_dia = TTAlgorithm(flow_asc_dia,freq,algorithm,waveform,continuous,show,handles.axes_asc_dia,flag_plot);
        transit_time_asc_dia = transit_time_asc_desc + transit_time_desc_dia;        
        TT_asc_dia_str = sprintf('%.5f', transit_time_asc_dia);
        set(handles.str_TT_asc_dia, 'String', TT_asc_dia_str);
        axes(handles.axes_asc_dia), ylim([ymin*1.2 ymax*1.2])
    case 'Least_squares_button'    
        arrayfun(@cla,findall(0,'type','axes')),set(handles.str_TT_asc_desc, 'String', []),set(handles.str_TT_desc_dia, 'String', []),set(handles.str_TT_asc_dia, 'String', []);
        algorithm = 3;
        %- Ascending-Descending
        flow_asc_desc = [flow(:,1) flow(:,2)];
        flag_plot = 'asc_desc';
        transit_time_asc_desc = TTAlgorithm(flow_asc_desc,freq,algorithm,waveform,continuous,show,handles.axes_asc_desc,flag_plot);
        TT_asc_desc_str = sprintf('%.5f', transit_time_asc_desc);
        set(handles.str_TT_asc_desc, 'String', TT_asc_desc_str);
        axes(handles.axes_asc_desc), ylim([ymin*1.2 ymax*1.2]), 
        %- Descending-Diaphragm
        flow_desc_dia = [flow(:,2) flow(:,3)];
        flag_plot = 'desc_dia';
        transit_time_desc_dia = TTAlgorithm(flow_desc_dia,freq,algorithm,waveform,continuous,show,handles.axes_desc_dia,flag_plot);
        TT_desc_dia_str = sprintf('%.5f', transit_time_desc_dia);
        set(handles.str_TT_desc_dia, 'String', TT_desc_dia_str);
        axes(handles.axes_desc_dia), ylim([ymin*1.2 ymax*1.2])
        %- Ascending-Diaphragm
        flow_asc_dia = [flow(:,1) flow(:,3)];
        flag_plot = 'asc_dia';
        transit_time_asc_dia = TTAlgorithm(flow_asc_dia,freq,algorithm,waveform,continuous,show,handles.axes_asc_dia,flag_plot);
        transit_time_asc_dia = transit_time_asc_desc + transit_time_desc_dia;
        TT_asc_dia_str = sprintf('%.5f', transit_time_asc_dia);
        set(handles.str_TT_asc_dia, 'String', TT_asc_dia_str);
        axes(handles.axes_asc_dia), ylim([ymin*1.2 ymax*1.2])
    case 'Cross_correlation_button'
        arrayfun(@cla,findall(0,'type','axes')),set(handles.str_TT_asc_desc, 'String', []),set(handles.str_TT_desc_dia, 'String', []),set(handles.str_TT_asc_dia, 'String', []);
        algorithm = 4;
        %- Ascending-Descending
        flow_asc_desc = [flow(:,1) flow(:,2)];
        flag_plot = 'asc_desc';
        transit_time_asc_desc = TTAlgorithm(flow_asc_desc,freq,algorithm,waveform,continuous,show,handles.axes_asc_desc,flag_plot);
        TT_asc_desc_str = sprintf('%.5f', transit_time_asc_desc);
        set(handles.str_TT_asc_desc, 'String', TT_asc_desc_str);
        axes(handles.axes_asc_desc), ylim([ymin*1.2 ymax*1.2]), 
        %- Descending-Diaphragm
        flow_desc_dia = [flow(:,2) flow(:,3)];
        flag_plot = 'desc_dia';
        transit_time_desc_dia = TTAlgorithm(flow_desc_dia,freq,algorithm,waveform,continuous,show,handles.axes_desc_dia,flag_plot);
        TT_desc_dia_str = sprintf('%.5f', transit_time_desc_dia);
        set(handles.str_TT_desc_dia, 'String', TT_desc_dia_str);
        axes(handles.axes_desc_dia), ylim([ymin*1.2 ymax*1.2])
        %- Ascending-Diaphragm
        flow_asc_dia = [flow(:,1) flow(:,3)];
        flag_plot = 'asc_dia';
        transit_time_asc_dia = TTAlgorithm(flow_asc_dia,freq,algorithm,waveform,continuous,show,handles.axes_asc_dia,flag_plot);
        transit_time_asc_dia = transit_time_asc_desc + transit_time_desc_dia;
        TT_asc_dia_str = sprintf('%.5f', transit_time_asc_dia);
        set(handles.str_TT_asc_dia, 'String', TT_asc_dia_str);
        axes(handles.axes_asc_dia), ylim([ymin*1.2 ymax*1.2])
end
%% Update flow analysis axis
Q = flow(:,1); t = [1:length(Q)]./freq;
[c i_systole] = max(Q);
Qt = Q(i_systole:end);
ind_neg = intersection(Qt);
if isempty(ind_neg)
    [c ind] = findpeaks(-Qt);
    ind_intersection = ind(1);    
else
    [ind_intersection] = ind_neg(1);
end
i_diastole = i_systole + ind_intersection(1) - 1;
Qes = Q(1:i_systole); tes = t(1:i_systole);
Qls = Q(i_systole:i_diastole);tls = t(i_systole:i_diastole);
%- Get tNd
t_pre_ejection = get_t_pre_ejection(Q);
tNd = t_pre_ejection/i_diastole;
%- Plot and EF1
axes(handles.axes_flow_analysis)
hold all,
hes = area(tes,Qes); hes.FaceColor = [.9 0 0]; hes.FaceAlpha = 0.2;
hls = area(tls,Qls); hls.FaceColor = [0 0 .9]; hls.FaceAlpha = 0.2;
plot(t,Q,'b','linewidth',2);
ylabel('Q (ml/s)'),xlabel('t (s)')
Area_es = trapz(tes,Qes);Area_sys = trapz(t(1:i_diastole),Q(1:i_diastole));
EF1 = Area_es*100/Area_sys;
EF1_str = sprintf('%.0f', EF1);
set(handles.str_EF1, 'String', EF1_str);
setappdata(0,'EF1',EF1);
%- HR
HR = 60/t(end);
HR_str = sprintf('%.1f', HR);
set(handles.str_HR, 'String', HR_str);
setappdata(0,'HR',HR);
%- SV
SV = trapz(t,Q);
SV_str = sprintf('%.0f', SV);
set(handles.str_SV, 'String', SV_str);
setappdata(0,'SV',SV);
%- CO
CO = SV*1e-3*HR;
CO_str = sprintf('%.1f', CO);
set(handles.str_CO, 'String', CO_str);
setappdata(0,'CO',CO);
%- Distensibility
if isempty(PP)==0;
    %- Ascending
    Amin_asc = min(areamm2(:,1)); Amax_asc = max(areamm2(:,1));
    Dis_Asc = (Amax_asc-Amin_asc)/Amin_asc/PP;
    Dis_Asc_str = sprintf('%.1f', Dis_Asc);
    set(handles.str_dis_asc, 'String', Dis_Asc_str);
    %- Descending
    Amin_desc = min(areamm2(:,2)); Amax_desc = max(areamm2(:,2));
    Dis_Desc = (Amax_desc-Amin_desc)/Amin_desc/PP;
    Dis_Desc_str = sprintf('%.1f', Dis_Desc);
    set(handles.str_dis_desc, 'String', Dis_Desc_str);
    %- Diaphram
    Amin_dia = min(areamm2(:,3)); Amax_dia = max(areamm2(:,3));
    Dis_Dia = (Amax_dia-Amin_dia)/Amin_dia/PP;
    Dis_Dia_str = sprintf('%.1f', Dis_Dia);
    set(handles.str_dis_dia, 'String', Dis_Dia_str);
    return
else
    Dis_Asc = [];
    Dis_Desc = [];
    Dis_Dia = [];
end
setappdata(0,'Dis_Asc',Dis_Asc);
setappdata(0,'Dis_Desc',Dis_Desc);
setappdata(0,'Dis_Dia',Dis_Dia);
% Update handles structure
guidata(hObject, handles);

function str_PP_Callback(hObject, eventdata, handles)
global Amin_asc Amax_asc Amin_desc Amax_desc Amin_dia Amax_dia areamm2...
    Dis_Asc Dis_Desc Dis_Dia EDV EF PP DBP SBP Pes ENd_est Ees ENd_avg tNd SV
str_PP = get(hObject,'String')
PP = str2double(str_PP);
setappdata(0,'PP',PP);
%- Prepare values
if isempty(SBP)==0
    DBP = SBP-PP;
    DBP_str = sprintf('%.0f',DBP);
    set(handles.str_DBP,'string',DBP_str);
    setappdata(0,'DBP',DBP);
end
if isempty(DBP)==0
    SBP = DBP+PP;
    SBP_str = sprintf('%.0f',SBP);
    set(handles.str_SBP,'string',SBP_str);
    setappdata(0,'SBP',SBP);
    if isempty(Pes)==1
       Pes = 0.9*SBP; 
       Pes_str = sprintf('%.0f',Pes);
       set(handles.str_Pes,'string',Pes_str);
       setappdata(0,'Pes',Pes);
    end
end

%- Compute Elastance at the onset of ejection (ENd_est)
if isempty(EF)==0 && isempty(DBP)==0 && isempty(Pes)==0
    EF = EF/100;
    ENd_avg = compute_ENd_avg(tNd);
    ENd_est = compute_ENd_est(EF,DBP,Pes,ENd_avg);
    ENd_est_str = sprintf('%.2f', ENd_est);
    set(handles.str_ENd_est,'string',ENd_est_str);
    setappdata(0,'ENd_est',ENd_est);
    EF = EF*100;
end
%- Compute end-systolic Elastance (Ees)
if isempty(DBP)==0 && isempty(ENd_est)==0 && isempty(SBP)==0 && isempty(SV)==0
   Ees = compute_Ees(DBP,ENd_est,SBP,SV);
   Ees_str = sprintf('%.2f', Ees);
   set(handles.str_Ees,'string',Ees_str);
   setappdata(0,'Ees',Ees);
end
%- Compute Distensibility 
if isempty(areamm2)==0
   [Dis_Asc Dis_Desc Dis_Dia] = compute_distensibility(PP,areamm2);
   Dis_Asc_str = sprintf('%1.1e', Dis_Asc); Dis_Desc_str = sprintf('%1.1e', Dis_Desc); Dis_Dia_str = sprintf('%1.1e', Dis_Dia);
   set(handles.str_dis_asc, 'String', Dis_Asc_str);set(handles.str_dis_desc, 'String', Dis_Desc_str);set(handles.str_dis_dia, 'String', Dis_Dia_str);
   setappdata(0,'Dis_Asc',Dis_Asc);setappdata(0,'Dis_Desc',Dis_Desc);setappdata(0,'Dis_Dia',Dis_Dia);
end

% --- Executes during object creation, after setting all properties.
function str_PP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function str_DBP_Callback(hObject, eventdata, handles)
global Amin_asc Amax_asc Amin_desc Amax_desc Amin_dia Amax_dia areamm2...
    Dis_Asc Dis_Desc Dis_Dia EDV EF PP DBP SBP Pes ENd_est Ees SV tNd ENd_avg
str_DBP = get(hObject,'String')
DBP = str2double(str_DBP);
setappdata(0,'DBP',DBP);
%- Prepare values
if isempty(SBP)==0
    PP = SBP-DBP;
    PP_str = sprintf('%.0f',PP);
    set(handles.str_PP,'string',PP_str);
    setappdata(0,'PP',PP);
end
if isempty(PP)==0
    SBP = DBP+PP;
    SBP_str = sprintf('%.0f',SBP);
    set(handles.str_SBP,'string',SBP_str);
    setappdata(0,'SBP',SBP);
    if isempty(Pes)==1
        Pes = 0.9*SBP;
        Pes_str = sprintf('%.0f',Pes);
        set(handles.str_Pes,'string',Pes_str);
        setappdata(0,'Pes',Pes);
    end
end
%- Compute Elastance at the onset of ejection (ENd_est)
if isempty(EF)==0 && isempty(DBP)==0 && isempty(Pes)==0
    EF = EF/100;
    ENd_avg = compute_ENd_avg(tNd);
    ENd_est = compute_ENd_est(EF,DBP,Pes,ENd_avg);
    ENd_est_str = sprintf('%.2f', ENd_est);
    set(handles.str_ENd_est,'string',ENd_est_str);
    setappdata(0,'ENd_est',ENd_est);
    EF = EF*100;
end
%- Compute end-systolic Elastance (Ees)
if isempty(DBP)==0 && isempty(ENd_est)==0 && isempty(SBP)==0 && isempty(SV)==0
   Ees = compute_Ees(DBP,ENd_est,SBP,SV);
   Ees_str = sprintf('%.2f', Ees);
   set(handles.str_Ees,'string',Ees_str);
   setappdata(0,'Ees',Ees);
end
%- Compute Distensibility 
if isempty(PP)==0 && isempty(areamm2)==0
   [Dis_Asc Dis_Desc Dis_Dia] = compute_distensibility(PP,areamm2);
   Dis_Asc_str = sprintf('%1.1e', Dis_Asc); Dis_Desc_str = sprintf('%1.1e', Dis_Desc); Dis_Dia_str = sprintf('%1.1e', Dis_Dia);
   set(handles.str_dis_asc, 'String', Dis_Asc_str);set(handles.str_dis_desc, 'String', Dis_Desc_str);set(handles.str_dis_dia, 'String', Dis_Dia_str);
   setappdata(0,'Dis_Asc',Dis_Asc);setappdata(0,'Dis_Desc',Dis_Desc);setappdata(0,'Dis_Dia',Dis_Dia);
end

% --- Executes during object creation, after setting all properties.
function str_DBP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function str_SBP_Callback(hObject, eventdata, handles)
global Amin_asc Amax_asc Amin_desc Amax_desc Amin_dia Amax_dia areamm2...
Dis_Asc Dis_Desc Dis_Dia EDV EF PP DBP SBP Pes ENd_est Ees SV tNd ENd_avg
str_SBP = get(hObject,'String')
SBP = str2double(str_SBP);
setappdata(0,'SBP',SBP);
%- Prepare values
if isempty(DBP)==0
    PP = SBP-DBP;
    PP_str = sprintf('%.0f',PP);
    set(handles.str_PP,'string',PP_str);
    setappdata(0,'PP',PP);
end
if isempty(Pes)==1
    Pes = SBP*0.9;
    Pes_str = sprintf('%.0f',Pes);
    set(handles.str_Pes,'string',Pes_str);
    setappdata(0,'Pes',Pes);
end
if isempty(PP)==0 
    DBP = SBP-PP;
    DBP_str = sprintf('%.0f',DBP);
    set(handles.str_DBP,'string',DBP_str);
    setappdata(0,'DBP',DBP);
end
%- Compute Elastance at the onset of ejection (ENd_est)
if isempty(EF)==0 && isempty(DBP)==0 && isempty(Pes)==0
    EF = EF/100;
    ENd_avg = compute_ENd_avg(tNd);
    ENd_est = compute_ENd_est(EF,DBP,Pes,ENd_avg);
    ENd_est_str = sprintf('%.2f', ENd_est);
    set(handles.str_ENd_est,'string',ENd_est_str);
    setappdata(0,'ENd_est',ENd_est);
    EF = EF*100;
end
%- Compute end-systolic Elastance (Ees)
if isempty(DBP)==0 && isempty(ENd_est)==0 && isempty(SBP)==0 && isempty(SV)==0
   Ees = compute_Ees(DBP,ENd_est,SBP,SV);
   Ees_str = sprintf('%.2f', Ees);
   set(handles.str_Ees,'string',Ees_str);
   setappdata(0,'Ees',Ees);
end
%- Compute Distensibility 
if isempty(PP)==0 && isempty(areamm2)==0
   [Dis_Asc Dis_Desc Dis_Dia] = compute_distensibility(PP,areamm2);
   Dis_Asc_str = sprintf('%1.1e', Dis_Asc); Dis_Desc_str = sprintf('%1.1e', Dis_Desc); Dis_Dia_str = sprintf('%1.1e', Dis_Dia);
   set(handles.str_dis_asc, 'String', Dis_Asc_str);set(handles.str_dis_desc, 'String', Dis_Desc_str);set(handles.str_dis_dia, 'String', Dis_Dia_str);
   setappdata(0,'Dis_Asc',Dis_Asc);setappdata(0,'Dis_Desc',Dis_Desc);setappdata(0,'Dis_Dia',Dis_Dia);
end

% --- Executes during object creation, after setting all properties.
function str_SBP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function str_Pes_Callback(hObject, eventdata, handles)
global Amin_asc Amax_asc Amin_desc Amax_desc Amin_dia Amax_dia areamm2...
Dis_Asc Dis_Desc Dis_Dia EDV EF PP DBP SBP Pes ENd_est Ees SV tNd ENd_avg
str_Pes = get(hObject,'String')
Pes = str2double(str_Pes);
setappdata(0,'Pes',Pes);
%- Compute Elastance at the onset of ejection (ENd_est)
if isempty(EF)==0 && isempty(DBP)==0 && isempty(Pes)==0
    EF = EF/100;
    ENd_avg = compute_ENd_avg(tNd);
    ENd_est = compute_ENd_est(EF,DBP,Pes,ENd_avg);
    ENd_est_str = sprintf('%.2f', ENd_est);
    set(handles.str_ENd_est,'string',ENd_est_str);
    setappdata(0,'ENd_est',ENd_est);
    EF = EF*100;
end
%- Compute end-systolic Elastance (Ees)
if isempty(DBP)==0 && isempty(ENd_est)==0 && isempty(SBP)==0 && isempty(SV)==0
   Ees = compute_Ees(DBP,ENd_est,SBP,SV);
   Ees_str = sprintf('%.2f', Ees);
   set(handles.str_Ees,'string',Ees_str);
   setappdata(0,'Ees',Ees);
end

% --- Executes during object creation, after setting all properties.
function str_Pes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function str_EF_Callback(hObject, eventdata, handles)
global Amin_asc Amax_asc Amin_desc Amax_desc Amin_dia Amax_dia areamm2...
Dis_Asc Dis_Desc Dis_Dia EDV EF PP DBP SBP Pes ENd_est Ees SV tNd ENd_avg
str_EF = get(hObject,'String')
EF = str2double(str_EF);
setappdata(0,'EF',EF);
%- Compute value
if isempty(SV)==0
    EF = EF/100;
    EDV = SV/EF;
    EDV_str = sprintf('%.0f',EDV);
    set(handles.str_EDV,'string',EDV_str);
    setappdata(0,'EDV',EDV);
    EF = EF*100;
end
%- Compute Elastance at the onset of ejection (ENd_est)
if isempty(EF)==0 && isempty(DBP)==0 && isempty(Pes)==0
    EF = EF/100;
    ENd_avg = compute_ENd_avg(tNd);
    ENd_est = compute_ENd_est(EF,DBP,Pes,ENd_avg);
    ENd_est_str = sprintf('%.2f', ENd_est);
    set(handles.str_ENd_est,'string',ENd_est_str);
    setappdata(0,'ENd_est',ENd_est);
    EF = EF*100;
end
%- Compute end-systolic Elastance (Ees)
if isempty(DBP)==0 && isempty(ENd_est)==0 && isempty(SBP)==0 && isempty(SV)==0
   Ees = compute_Ees(DBP,ENd_est,SBP,SV);
   Ees_str = sprintf('%.2f', Ees);
   set(handles.str_Ees,'string',Ees_str);
   setappdata(0,'Ees',Ees);
end

% --- Executes during object creation, after setting all properties.
function str_EF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function str_EDV_Callback(hObject, eventdata, handles)
global Amin_asc Amax_asc Amin_desc Amax_desc Amin_dia Amax_dia areamm2...
Dis_Asc Dis_Desc Dis_Dia EDV EF PP DBP SBP Pes ENd_est Ees SV tNd ENd_avg
str_EDV = get(hObject,'String')
EDV = str2double(str_EDV);
setappdata(0,'EDV',EDV);
%- Create value
if isempty(SV)==0
    EF = round(SV*100/EDV);
    EF_str = sprintf('%.0f',EF);
    set(handles.str_EF,'string',EF_str);
    setappdata(0,'EF',EF);
end
%- Compute Elastance at the onset of ejection (ENd_est)
if isempty(EF)==0 && isempty(DBP)==0 && isempty(Pes)==0
    EF = EF/100;
    ENd_avg = compute_ENd_avg(tNd);
    ENd_est = compute_ENd_est(EF,DBP,Pes,ENd_avg);
    ENd_est_str = sprintf('%.2f', ENd_est);
    set(handles.str_ENd_est,'string',ENd_est_str);
    setappdata(0,'ENd_est',ENd_est);
    EF = EF*100;
end
%- Compute end-systolic Elastance (Ees)
if isempty(DBP)==0 && isempty(ENd_est)==0 && isempty(SBP)==0 && isempty(SV)==0
   Ees = compute_Ees(DBP,ENd_est,SBP,SV);
   Ees_str = sprintf('%.2f', Ees);
   set(handles.str_Ees,'string',Ees_str);
   setappdata(0,'Ees',Ees);
end

% --- Executes during object creation, after setting all properties.
function str_EDV_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Continue_button.
function Continue_button_Callback(hObject, eventdata, handles)
global transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia
setappdata(0,'transit_time_asc_desc',transit_time_asc_desc)
setappdata(0,'transit_time_desc_dia',transit_time_desc_dia)
setappdata(0,'transit_time_asc_dia',transit_time_asc_dia)
uiresume


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADDITIONAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TT = TTAlgorithm(signal,f,algorithm,waveform,continuous,show,handles_axes,flag_plot)
% Description:
%   TTAlgorithm applies user defined pulse wave analysis algorithms in
%   order to determine the transit time between two wave forms.  This
%   software assumes that the two waveforms are measured simultaneously,
%   and no temporal alignment is needed.setappdata(0,'time_transit',transit_time)
% Features:
% - Foot to foot algorithm
% - Foot to foot radius algorithm
% - Cross correlation of the entire cycle
% - Least Squares
% 
% signal - matrix with two physiological waveform vectors
% f - signal freqency (Hz), (note the frequency of the Test Doppler Data was 100Hz)
% algorithm - 1, 2, 3, or 4 refering to the above algorithms respectively
% waveform - 1, or 2, referring to pressure/area or velocity/flow respectively
% continuous - 0, or 1, indicating a single wave (0) or multiple waves (1)
% show - display location of 'feet' used for analysis, (1=yes, 0=no)
% TT - vector of transit times, length depends on number of waveforms
% 
% ************De
% ** Please cite the following paper if used for publication **
% 
% N R Gaddum, J Alastruey, P Beerbaum, P Chowienczyk, T Schaeffter,
% A technical assessment of pulse wave velocity algorithms applied 
% to non-invasive arterial waveforms, Annals of Biomedical Engineering,
% 2013, 41(12):2617-29. DOI: 10.1007/s10439-013-0854-ysetappdata(0,'time_transit',transit_time)
% 
% ************
% 
% Author:
% Dr. Nicholas Gaddum
% 
% Revisions:
% 2012-Aug-06   Created function
% 2012-Aug-09   Removed replication of time parameters.  Only frequency 
% 	remaPWV_Calulator_Least_Squaresining.  Added pressure/area vs. velocity/flow parameter for
% 	Find_Landmarks observer.  Bug fixing.
% 2012-Aug-13   Revised README file.
% 2012-Aug-21   Pressure test data added to the Test Data File.  Added 
% 	recomendations.  Changes to the scanning window for algorithm 1.
% 2013-Aug-10   Bug fixing for high resolution pressure/flow data. Higher
% 	stability for foot location using the minimum radius
% 2014-Mar-24   Updated README file.  Added troubleshooting section
% 2014-Oct-02   Reduction of m files, and include option for single or
%   multiple waveform analysis
% 2014-Dec-03   Bug fix for brachial-ankle, carotid-ankle data, where the
% 	feet have a greater spacing.  Cross correlation bug fix
% 2014-Dec-19 	Removed and unused .mat file
% 2016-Jun-10   Bug fixing for low temporal resolution MRI flow data, and
% 	added signal observers for better stability of Cross Correlation

if nargin < 1
    error('Please specify a matrix of two waveform vectors')
    return
elseif nargin == 1
    warning('Sample frequency, algorithm, and data type not chosen.  Set to defaults, (f=100Hz, Foot to foot, pressure data)')
    f = 100;
    algorithm = 1;
    waveform = 1;
    show = 0;
elseif nargin == 2
    warning('Algorithm and data type not chosen.  Set to defaults, (Foot to foot, pressure data)')
    algorithm = 1;
    waveform = 1;
    show = 0;
elseif nargin == 3
    warning('Data type not chosen.  Set to default, (pressure data)')
    waveform = 1;
    show = 0;
elseif nargin == 4
    show = 0;
else
end

if algorithm > 4
    warning('Invalid "algorithm" code.  Set to default, (foot to foot)')
    algorithm = 1;
end
if waveform > 2
    warning('Invalid "waveform" code.  Set to default, (velocity)')
    waveform = 2;
end
if continuous > 1
    warning('Invalid "continuous" code.  Set to default, (signal wave)')
    continuous = 0;
end
if show > 1
    warning('Invalid "show" code.  Set to default, (show plots)')
    show = 1;
end
if f < 50
    warning('Check prescribed frequency, it appears to be too low')
end

if show == 1
    if algorithm == 1
        fprintf('Using Foot to foot:\n')
    elseif algorithm == 2
        fprintf('Using Foot to foot minimum radius:\n')
    elseif algorithm == 3
        fprintf('Using Least Squares:\n')
    elseif algorithm == 4
        fprintf('Using Cross Correlation of the entire cycle:\n')
    end
end

% Prepare data for processing, ensure signals are in two rows
if size(signal,1) > size(signal,2)
    signal = signal';
end
ind = find(abs(signal) < 10^(-10));
signal(ind) = 10^(-10); % Replace zero data with near zero data, (for Cross correlation algorithm)

dt = 1/f;
n = size(signal,2);
t = [0:dt:(n-1)*dt];

%%%%%%% NOTE, wave feet can be poorly locaFind_Landmarksted if your frequency is
%%%%%%% significantly low/high.  In this case, try varying the value of the
%%%%%%% following two variables
npoly_1 = ceil(0.03*f); % number of terms taken in the linear polyfit
npoly1_1 = ceil(0.10*f); % number of terms taken in the trend radius calculations

if continuous == 1
    t_intv = Cycle_Approximator(t,signal);
elseif continuous == 0
    t_intv = 0;
    
    % If a single cycle is used, (particularly important for MRI), add some
    % extra data before the foot of each profile in order for the software
    % to be able to scan the waveform through the foot of the waveform
    n1 = 10;
    Amp1 = max(signal(1,:)) - min(signal(1,:));
    m1 = Amp1 * 0.001;
    Amp2 = max(signal(2,:)) - min(signal(2,:));
    m2 = Amp2 * 0.001;
    noise = rand(1,n1) - 0.5;
    trash = [signal(1,1)*ones(1,n1); signal(2,1)*ones(1,n1)] + [m1*noise; m2*noise];
    
    signal = [trash, signal];
    t = [0:dt:(n-1+n1)*dt];
end

% Find maxima/minima/maximum gradients for all algorithms
[indmaxs,indmins,gradins] = Find_Landmarks(t,signal,t_intv,npoly_1,npoly1_1,f,algorithm,waveform,continuous,show,handles_axes,flag_plot);

switch algorithm
    case 3
        % Least Squares difference of the systolic upstroke
        t_int = PWV_Calulator_Least_Squares(t,signal,indmaxs,indmins,handles_axes);
    case 4
        % Cross correlation of the entire waveform
        t_int = PWV_Calulator_CC_cycle(t,signal,indmaxs,indmins,gradins,show,handles_axes);
    otherwise
        % Both foot to foot algorithms
        t_int = PWV_Calulator_FTF(t,signal,indmaxs,indmins,gradins,show,handles_axes);
end
TT = t_int;

function [t_int,t_foot,horint] = PWV_Calulator_FTF(t,signal,indmaxs,indmins,gradins,show,handles_axes)

%%%%%%%%%%%%% Find the 'Feet - Start %%%%%%%%%%%%%
lin = min(length(indmaxs),size(indmins,2));
lin = min(lin,size(gradins,2));

for i=1:size(signal,1)
    for j = 1:size(gradins,2)
        if gradins(i,j,1) == 0 || indmins(i,j,1) == 0
        else
            % Horizontal limits
            centre = indmins(i,j,1);
            space = 1;
            horint(i,j) = mean( signal(i, (centre-space):(centre+space) ) ); % Add one to align minima and gradient coordinates
            
            % Locate trasnient intercept
            t_foot(i,j) = (horint(i,j) - gradins(i,j,3))/(gradins(i,j,2));
            ind(i,j) = gradins(i,j,1);
        end
    end
end
TT = t_foot(2,:) - t_foot(1,:);
if median(TT) < 0
    TT = -TT;
end

switch show
    case 0
    case 1
        axes(handles_axes)
        plot(t_foot,horint,'rx','MarkerSize',12,'Linewidth',2)
end

non_zeros = find(t_foot(1,:).*t_foot(2,:));

t_int = TT;

function t_int= PWV_Calulator_Least_Squares(t,signal,indmaxs,indmins,handles_axes)

signal = signal;
t = t;

%%%%%%%%%%%%% Least Squares - Start %%%%%%%%%%%%%

dt = t(2)-t(1);

trash = diff(indmaxs);
trash = median(trash);

ind = find(indmins(1,:,1));
kernel1_tot = [indmins(1, : ,1);indmaxs(1, : )];
ind = find(indmins(2,:,1));
kernel2_tot = [indmins(2, : ,1);indmaxs(2, : )];
Prox = 4;

pixPCycle = round(abs(mean(diff(kernel1_tot(1,:)))));

padding1 = round(pixPCycle/3);
padding2 = round(2*pixPCycle/3);

count = 1;

if length(ind) > 1 % if there are more than one cycles
    
    for k = 2:length(ind)
        
        kernel1 = kernel1_tot(:,k);
        kernel2 = kernel2_tot(:,k);
        
        if kernel1(1)<kernel1(2) && kernel2(1)<kernel2(2)
            if kernel1(1) > kernel2(1)
                kernel1(1) = kernel2(1);
            end
            shift = round( abs(kernel1(2) - kernel1(1)));
            
            max_lengthBack = max( [(kernel1(2) - kernel1(1)) , (kernel2(2) - kernel2(1)) , abs(kernel2(2) - kernel1(1))] );
            max_lengthForward = round(max_lengthBack/2);
            
            if size(signal,2) < kernel2(2)+max_lengthForward
                max_lengthForward = round(size(signal,2) - kernel2(2));
            end
            
            % Normalise the signals and lift off the x axis
            min1 = min(signal(1,kernel1(1):kernel1(2)));
            trash1 = signal(1,:) - min1;
            sig1_norm = trash1/(trash1(kernel1(2))) + 0.05;
            
            min2 = min(signal(2,kernel2(1):kernel2(2)));
            trash1 = signal(2,:) - min2;
            sig2_norm = trash1/(trash1(kernel2(2)));
            
            sig2_base = zeros(1,length(t));
            lim1 = kernel2(2)-max_lengthBack;
            if lim1 < 1
                lim1 = 1;
            end
            lim2 = kernel2(2)+max_lengthForward;
            if lim2 > length(sig2_base)
                lim2 = length(sig2_base);
            end
            sig2_base(lim1:lim2) = sig2_norm(lim1:lim2) + 0.05;
            
            if (kernel1(1) - padding1) < 1
                padding1 = kernel1(1) - 1;
            elseif (kernel1(2) + padding2) >= length(t)
                padding2 = length(t) - kernel1(2);
            end
            
            % Signals for least squares correlation
            sig1_new = sig1_norm(kernel1(1)-padding1 : kernel1(2)+padding2);
            sig2     = sig2_base(kernel1(1)-padding1 : kernel1(2)+padding2);
            
            CCLength = length(find(sig2));
            
            y = [];
            x = [];
            
            loop = 1;
            j = -round(shift);
            
            while loop == 1
                if j<0
                    sig2_new = [zeros(1,-j),sig2(1 : end-(-j)) ];
                else
                    sig2_new = [sig2(j+1:end),zeros(1,j)];
                end
                
                trash = find(sig2_new);
                if (trash(end)-trash(1)) >= round(CCLength*0.7)
                    trash = (sig1_new.*sig2_new);
                    [C,I] = find(trash);
                    LS = 0;
                    if numel(I)>0
                        LS = sum( (sig1_new(I) -sig2_new(I)).^2 );
                    else
                        LS = 20;
                    end
                    
                    y = [y,LS];
                    x = [x,(j*dt)];
                    
                    j = j+1;
                    ind = find(sig2_new);
                    if ind < (max_lengthBack+1)
                        loop = 0;
                    end
                else
                    if isempty(y)
                        j = j+1;
                    else
                        loop = 0;
                        skip = 1;
                    end
                end
            end
            
            intZero = find(x==0);
            
            if ~isempty(intZero)
                
                trash1 = diff(y);
                MaxMinInd = find(trash1(1:end-1).*trash1(2:end)<0)+1;
                if length(MaxMinInd)>1
                    MaxMin = y(MaxMinInd);
                    startGrad = diff(MaxMin(1:2));
                    if startGrad > 0
                        trash2 = MaxMinInd(1:2:end);
                    else
                        trash2 = MaxMinInd(2:2:end);
                    end
                else
                    trash2 = MaxMinInd;
                end
                
                trash3 = abs( (intZero+Prox - trash2));
                trash4 = (trash3) .* y(trash2);
                [~,ind] = min(trash4);
                I = trash2(ind);
                
                interp_span = 1;
                if (I+interp_span) > length(y)
                    I = length(y) - interp_span;
                elseif (I-interp_span) < 1
                    I = 1 + interp_span;
                end
                interp_int = [I-interp_span:I+interp_span];
                int = find(interp_int>0);
                interp_int = interp_int(int);
                
                poly = polyfit(x(interp_int),y(interp_int),2);
                x1 = [x(1):0.001:x(end)];
                y1 = poly(1).*x1.^2 + poly(2).*x1 + poly(3);
                
                m_TT = -poly(2)/(2*poly(1));
                
                t_int(count) = m_TT;
                count = count+1;
                skip = 0;
            else
                t_int(count) = nan;
                count = count+1;
            end
        end
    end
    setappdata(0,'time_transit',transit_time)
else
    
    ind = indmins(:,:,1);
    kernel1 = [indmins(1,:,1) , indmaxs(1)];
    kernel2 = [indmins(2,:,1) , indmaxs(2)];
    
    shift = max( round( abs(kernel1(2) - kernel1(2))*5 ),20 );
    
    max_length = max( [(kernel1(2) - kernel1(1)) , (kernel2(2) - kernel2(1))] );
    
    trash1 = signal(1,:) - signal(1,kernel1(1));
    sig1_norm = trash1/(trash1(kernel1(2))) + 0.05;
    
    trash1 = signal(2,:) - signal(2,kernel2(1));
    sig2_norm = trash1/(trash1(kernel2(2)));
    
    sig2_base = zeros(1,length(t));
    sig2_base(kernel2(2)-max_length:kernel2(2)) = sig2_norm(kernel2(2)-max_length:kernel2(2)) + 0.05;
    
    y = [];
    x = [];
    
    loop = 1;
    j = -round(shift/4);
    while loop == 1
        if j<0
            sig2_new = [zeros(1,-j),sig2_base(1 : end-(-j)) ];
        else
            sig2_new = [sig2_base(j+1:end),zeros(1,j)];
        end
        
        trash = find(sig2_new);
        if numel(trash)<1
            loop = 0;
        elseif (trash(end)-trash(1)) >= (kernel2(2) - kernel2(1)-1)
            trash = (sig1_norm.*sig2_new);
            [C,I] = find(trash);
            LS = 0;
            if numel(I)>0
                LS = sum( (sig1_norm(I) -sig2_new(I)).^2 );
            else
                LS = 20;
            end
            
            y = [y,LS];
            x = [x,(j*dt)];
            
            j = j+1;
            ind = find(sig2_new);
            if ind < (max_length+1)
                loop = 0;
            end
        else
        end
    end
    
    [~,I] = min(y);
    interp_span = 5;
    if (I+interp_span) > length(y)
        I = length(y) - interp_span;
    elseif (I-interp_span) < 1
        I = 1 + interp_span;
    end
    interp_int = [I-interp_span:I+interp_span];
    int = find(interp_int>0);
    interp_int = interp_int(int);
    
    poly = polyfit(x(interp_int),y(interp_int),2);
    
    t_int = -poly(2)/(2*poly(1));
end

function t_int = PWV_Calulator_CC_cycle(t,signal,indmaxs,indmins,gradins,show,handles_axes)
signal = signal;
t = t;

dt = t(2)-t(1);

shift = round(0.2/dt);

TTApprox = median(diff(gradins(:,:,1)));

axes(handles_axes)
hold on

if size(gradins,2) > 1
    cla
    
    t_int = [];
    
    kernel = gradins(:,:,1);
    if mean(kernel(1,:) - kernel(2,:))<0
    else
        signal = flipud(signal);
    end
    
    % Determine which of the two velocity profiles have the most indexed max
    % gradients
    trash_1 = find(kernel(1,:));
    trash_2 = find(kernel(2,:));
    if length(trash_1) >= length(trash_2)
        kernel = kernel(1,:);
    else
        kernel = kernel(2,:);
    end
    
    ind = find(kernel);
    
    for i=2:length(ind)
        
        
        kernel1 = [kernel(ind(i)) , kernel(ind(i-1))];
        
        if kernel1(1)-(2*shift) < 1
            shift = floor( (kernel1(1))/2 - 1 );
        end
        
        sig1 = signal(1,kernel1(1)-(2*shift) : kernel1(2)-shift);
        sig2 = signal(2,kernel1(1)-(2*shift) : kernel1(2)-shift);
        correlationEnvelope = length(sig1) - 2*TTApprox;
        
        % Normalise Carotid Wave Data
        trash = sig1 - min(sig1);
        sig1 = trash/max(trash);
        
        % Normalise Femoral Wave Data
        trash = sig2 - min(sig2);
        sig2 = trash/max(trash);
        
        y = [];
        x = [];
        for j= -round(TTApprox/4): 1 : 1.5*TTApprox
            
            sig1_new = sig1( (shift+1)   : (shift+correlationEnvelope)   );
            sig2_new = sig2( (shift+1)+j : (shift+correlationEnvelope)+j );
            
            subplot(1,2,1)
            cla
            hold on
            plot(sig1_new,'r')
            plot(sig2_new,'k')
            
            trash = find(sig2_new);
            
            xx = sum(sig1_new.*sig2_new);
            y = [y,xx];
            x = [x,(j*dt)];
            
            subplot(1,2,2)
            plot(x,y,'xk','MarkerSize',12)
            
            if xx > max(y(1:end-1))
                max_ind = j;
                max_CC = xx;
            else
            end
        end
        
        loop=1;
        k=1;
        while loop==1
            k=k+1;
            poly = polyfit(x(k-1:k+1),y(k-1:k+1),1);
            if poly(1) > 10
                ind1 = k;
                loop = 0;
            end
            if (k+1) == length(y)
                ind1 = 1;
                loop = 0;
            end
        end
        
        [~,I] = max(y(ind1:end));
        I = I + (ind1-1);
        interp_span = 5;
        
        if length(y) < (I+interp_span)
            interp_span = length(y) - I;
        elseif (I-interp_span) < 1
            interp_span = I-1;
        end
        
        interp_int = [I-interp_span:I+interp_span];
        int = find(interp_int);
        interp_int = interp_int(int);
        
        if var(y)<0.1
            nick = 1;
        end
        
        poly = polyfit(x(interp_int),y(interp_int),2);
        x1 = [x(1):0.001:x(end)];
        y1 = poly(1).*x1.^2 + poly(2).*x1 + poly(3);
        
        m_TT = -poly(2)/(2*poly(1));
        
        t_int = [t_int, m_TT];
    end
    
    
else
    
    
    ind = indmins(:,:,1);
    kernel1 = [indmins(1,:,1) , indmaxs(1)];
    kernel2 = [indmins(2,:,1) , indmaxs(2)];
    
    shift = max( round( abs(kernel2(2) - kernel1(2))*4 ) , 25 );
    
    trash = round(size(signal,2)/4);
    [~,I] = min( signal(1, (kernel1(1)+trash) :end));
    max_length = I  + kernel1(1)+trash - 1;
    
    sig1_new = signal(1,:);
    
    sig2 = zeros(1,size(signal,2));
    sig2(kernel2(1):max_length) = signal(2,kernel2(1):max_length);
    
    y = [];
    x = [];
    
    for j= -round(shift):1:round(shift)
        if j<0
            sig2_new = [zeros(1,-j),sig2(1 : end-(-j)) ];
        else
            sig2_new = [sig2(j+1:end),zeros(1,j)];
        end
        
        trash = find(sig2_new);
        if numel(trash)<1
            save('error_data','')
        elseif (trash(end)-trash(1)) >= (length( (kernel2(1)+1) :max_length)-1)
            xx = sum(sig1_new.*sig2_new);
            y = [y,xx];
            x = [x,(j*dt)];
            
            if xx > max(y(1:end-1))
                max_ind = j;
                max_CC = xx;
            else
            end
        else
        end
    end
    
    % Recreate max correlation to find the Correlation coefficient
    j = max_ind;
    if j<0
        sig2_new = [zeros(1,-j),sig2(1 : end-(-j)) ];
    else
        sig2_new = [sig2(j+1:end),zeros(1,j)];
    end
    
    ind_x = find(sig2_new);
    SumProx1 = sum(sig1_new(ind_x).^2);
    SumProx2 = sum(sig2_new(ind_x).^2);
    CC_ref = max(SumProx1,SumProx2);
    
    CrossCoeff = max_CC/CC_ref;
    
    loop=1;
    i=1;
    while loop==1
        i=i+1;
        poly = polyfit(x(i-1:i+1),y(i-1:i+1),1);
        if poly(1) > 10
            ind = i;
            loop = 0;
        end
        if (i+1) == length(y)
            ind = 1;
            loop = 0;
        end
    end
    
    [~,I] = max(y(ind:end));
    I = I + (ind-1);
    interp_span = 5;
    
    if length(y) < (I+interp_span)
        interp_span = length(y) - I;
    elseif (I-interp_span) < 1
        interp_span = I-1;
    end
    
    interp_int = [I-interp_span:I+interp_span];
    int = find(interp_int);
    interp_int = interp_int(int);
    
    poly = polyfit(x(interp_int),y(interp_int),2);
    x1 = [x(1):0.001:x(end)];
    y1 = poly(1).*x1.^2 + poly(2).*x1 + poly(3);
    
    t_int = -poly(2)/(2*poly(1));
    
%     switch show
%         case 0
%         case 1
%             plot(t_foot,horint,'xr','MarkerSize',12)
%     end
    
end

function [indmaxs,indmins,gradins] = Find_Landmarks(t,signal,t_intv,npoly_1,npoly1_1,f,algorithm,waveform,continuous,show,handles_axes,flag_plot)
%%%%%%%%%%%%% Locate Minima, Maxima and Max Gradients %%%%%%%%%%%%%
t = t';
signal = signal';

switch show
    case 0
    case 1
        axes(handles_axes)
        cla
        hold on
        if strcmp(flag_plot,'asc_desc')
            p1 = plot(t,signal(:,1),'b','linewidth',1);
            p2 = plot(t,signal(:,2),'r','linewidth',1);
            legend([p1 p2],{'Ascending','Descending'},'AutoUpdate','off')
            legend boxoff
        elseif strcmp(flag_plot,'desc_dia')
            p1 = plot(t,signal(:,1),'r','linewidth',1);
            p2 = plot(t,signal(:,2),'g','linewidth',1);
            legend([p1 p2],{'Descending','Diaphragm'},'AutoUpdate','off')
            legend boxoff
        else
            p1 = plot(t,signal(:,1),'b','linewidth',1);
            p2 = plot(t,signal(:,2),'g','linewidth',1);
            legend([p1 p2],{'Ascending','Diaphragm'},'AutoUpdate','off')
            legend boxoff
        end
%         low = min(min(signal));
%         high = max(max(signal));
        xlim([ t(1) t(end)]); %(low - 0.1*abs(low)) (high + 0.1*abs(high)) ])
        xlabel('Time (s)','FontSize',10)
        ylabel('Flow (ml/s)','FontSize',10)
        set(gca, 'Box', 'off' );
end

%%%%%%%%%%%%% Determine Minima, Maxima and Max Gradients - Start %%%%%%%%%%%%%%%%%
signal1 = signal';
t = t';

clear signal

for j = 1:size(signal1,1)
    %%%%% The Identification of Peaks %%%%%
    % indmax are the maxima positions, indmin are the minima positions
    indmax = nan;
    indmin = [];
    gradin = []; % Gradient info to store the index, gradient and y-intercept for the Find the Foot p-processing
    intercept = []; % plotted locations upon the polynomial of max gradient
    
    signal = signal1(j,:);
    
    if length(t_intv) > 1 % i.e. numerous cycles, (continuous = 1)
        threshold = 2.5.*sqrt(var( signal(t_intv(2) : t_int,v(1) ) ));
        if waveform == 1
            t_max_min = 0.50 * abs(mean(diff(t_intv)))/f;
        elseif waveform == 2
            t_max_min = 0.22 * abs(mean(diff(t_intv)))/f;
        end
    else % i.e. single cycle, (continuous = 0)
        threshold = 2.5.*sqrt(var( signal ));
        if waveform == 1
            t_max_min = 0.50 * length(signal)/f;
        elseif waveform == 2
            t_max_min = 0.22 * length(signal)/f;
        end
    end
    
    % Threshold should be set to the height of the smallest peak that should be
    % detected
    sa = min(signal);  % sa is the maximal element since the last minimum
    sb = max(signal);  % sb is the gradient since the last maximum
    sc = 0; % Start gradient off at a very high negative gradient, (i.e. mid/late systole)
    
    
    trash = diff(signal);
    equalAxisConst = (max(trash) / (t(2)-t(1))) / 10;
    
    i=length(signal)-2;
    a=length(signal);
    % d is the direction of the signal, 0-indetermined direction,
    % 1-increasing, 2-decreasing, (note this is observing the waveform from
    % the end to the beginning, i.e. t=end to t=t(1).
    if continuous == 0
        d=1; % this case assumes that only one waveform is available, i.e. MRI
    else
        d=0; % thsi case assumes continuous data with numerous waveforms
    end
    e=0;
    
    if t_intv(1) >= i
        count_thresh = 2;
    else
        count_thresh = 1;
    end
    
    v_b = [];
    
    while (i > npoly_1)
        npoly = npoly_1;
        npoly1 = npoly_1;
        i=i-1;
        
        if i == t_intv(count_thresh)%
            if (count_thresh+1) > length(t_intv)
                if t_intv(count_thresh+1)<1
                    count_thresh = count_thresh - 2;
                else
                end
                threshold = 2.5.*sqrt(var( signal(t_intv( count_thresh+1) : t_intv(count_thresh) ) ));
                count_thresh = count_thresh + 1;
            end
        else
        end
        
        if (d==0) % direction not yet determined
            if (sa >= (signal(i)+threshold))
                % d=2; % FALL
            elseif (signal(i) >= (sb+threshold))
                d=1; % RISE
            end;
            if (sa <= signal(i))
                sa = signal(i);
                a = i;
            elseif (signal(i) <= sb)
                sb = signal(i);
            end;
            
        elseif (d==1) % direction from min to max (RISE)
            if (sa <= signal(i))
                sa = signal(i);
                a = i;
            elseif (sa >= (signal(i)+threshold))% && ~isempty(indmin)
                % In order to prevent a repeating location of maxima
                % due to the regression of i to a, check for repitition
                % of maxima.
                if isnan(indmax) % I had to use an NAN starter, because I cannot get MATLAB to use isempty here
                    indmax = a;
                    i=a;
                    sb = inf;
                else
                    if a == indmax(end)
                    else
                        indmax = [indmax a];
                        i=a;
                        sb = inf;
                    end
                end
                sc = -10;
                b = i;
                d=2;
                e=2;
            end;
            
        elseif (d==2) && i < length(signal) || i < 6 %-limits1 && i > limits1 % direction from max to min (FALL)
            % Adjust window size for gradient and radius of curvature
            % approx.s near beginning/end of signal
            if i < (round( (npoly1+1)/2))
                npoly1 = 2*i-1;
            else
            end
            if i < (round( (npoly+1)/2))
                npoly = 2*i-1;
            else
            end
            
            if (i+ round( (npoly1+1)/2 )) > length(signal)
                npoly1 = 2 * (length(signal)-i)-1;
            else
            end
            if (i+ round( (npoly+1)/2 )) > length(signal)
                npoly = 2 * (length(signal)-i)-1;
            else
            end
            
            % To find the first derivative - for point of max gradient
            poly1 = polyfit(t(i+(round( (npoly+1)/2 )-1) :-1: i-(round( (npoly+1)/2 )-1)),...
                signal(i+(round( (npoly+1)/2 )-1) :-1: i-(round( (npoly+1)/2 )-1)),1);
            
            if algorithm > 1 % If definition of minimum is at the point of maximum radius of curvature
                % Use three points to estimate a radius of curvature
                y_1 = mean(signal(i+1 : i+(floor(npoly1/2))));
                x_1 = mean(     t(i+1 : i+(floor(npoly1/2)))) * equalAxisConst;
                y_2 = signal(i);
                x_2 = t(i) * equalAxisConst;
                y_3 = mean(signal(i-1 : -1 : i-(floor(npoly1/2))));
                x_3 = mean(     t(i-1 : -1 : i-(floor(npoly1/2)))) * equalAxisConst;
                
                % Form linear equations
                ma = (y_2-y_1) / (x_2-x_1);
                mb = (y_3-y_2) / (x_3-x_2);
                
                % Standardise the chord lengths
                x_5 = linspace(x_2,x_1,20);
                y_5 = ma.*(x_5 - x_2) + y_2;
                x_6 = linspace(x_2,x_3,20);
                y_6 = mb.*(x_6 - x_2) + y_2;
                l_5 = sqrt((x_5 - x_2).^2 + (y_5 - y_2).^2);
                l_6 = sqrt((x_6 - x_2).^2 + (y_6 - y_2).^2);
                diff_1 = min(l_5(end),l_6(end));
                % Redefine (x_1,y_1) and (x_3,y_3), (points either side of
                % signal(i)
                ind = find(l_5 > diff_1,1,'first');
                if isempty(ind)
                    ind = length(l_5);
                end
                x_1 = x_5(ind);
                y_1 = y_5(ind);
                ind = find(l_6 > diff_1,1,'first');
                if isempty(ind)
                    ind = length(l_6);
                end
                x_3 = x_6(ind);
                y_3 = y_6(ind);
                
                % Form linear equations
                ma = (y_2-y_1) / (x_2-x_1);
                mb = (y_3-y_2) / (x_3-x_2);
                a1 = [ (x_1+x_2) / 2 , (y_1+y_2) / 2 ];
                b1 = [ (x_2+x_3) / 2 , (y_2+y_3) / 2 ];
                
                % Coordinate of the centre of the LV arc
                x_4 = ( b1(2)-a1(2) - (a1(1)/ma - b1(1)/mb) ) / (-1/ma + 1/mb);
                y_4 = -1/ma * (x_4 - a1(1)) + a1(2);
                
                grad2 = ( ( (y_1-y_4)^2 + (x_1-x_4)^2 )^0.5 +...
                    ( (y_2-y_4)^2 + (x_2-x_4)^2 )^0.5 +...
                    ( (y_3-y_4)^2 + (x_3-x_4)^2 )^0.5 )/3;
                
                if (poly1(1) >= sc) && d==2
                    sc = poly1(1);
                    grad_store = [poly1(1),poly1(2)];
                    v = grad_store(1) * t(i) + grad_store(2);
                    c = [i ; poly1(1) ; poly1(2) ; v];
                else
                end
                
                % Check if a new minimum radius of curvature has been reached
                if (grad2 > 10^5*sb) || i < (npoly_1+1) || ((t(a) - t(i))>t_max_min) %&& && b~=size(signal,2))
                    min_data = [b ; v_b];
                    indmin = [indmin min_data];
                    sa = signal(i);
                    d=1;
                    
                    if e==2
                        gradin = [gradin c];
                        intercept = [intercept (c(2)*t(c(1)) + c(3))];
                        sc = 0;
                        e=1;
                        pause_grad = 1;
                    end
                elseif (grad2 <= sb) && y_4 > y_2 && i~=5 && (signal(a)-signal(i))>threshold/2 %(t(a) - t(i))>0.05
                    sb = grad2;
                    b = i;
                    % Mean velocity at t(i)
                    v_min = mean( signal(b+(round( (1+1)/2 )-1) :-1: b-(round( (1+1)/2 )-1)) );
                    v_b = v_min;
                end
                
            else % If definition of minimum is at the local minimum
                if (poly1(1) >= sc) && d==2
                    sc = poly1(1);
                    v = poly1(1) * t(i) + poly1(2);
                    c = [i ; poly1(1) ; poly1(2) ; v];
                else
                end
                
                % Check if a new minimum has been reached
                if signal(i) <= sb && ((t(a) - t(i))<t_max_min)
                    sb = signal(i);
                    b = i;
                    % Mean velocity at t(i)
                    v_min = mean( signal(b+1 :-1: b-1) );
                    v_b = v_min;
                elseif sb <= (signal(i)-threshold) || i < (npoly_1+1) || ((t(a) - t(i))>t_max_min)
                    min_data = [b ; v_b];
                    indmin = [indmin min_data];
                    sa = signal(i);
                    d=1;
                    
                    if e==2
                        gradin = [gradin c];
                        intercept = [intercept (c(2)*t(c(1)) + c(3))];
                        sc = 0;
                        e=1;
                        pause_grad = 1;
                    end
                end
            end
        end
    end
    
    if isnan(indmax)
        error('Could not determine features of the physiological waveform - please chack signal data...')
    end
    
    switch show
        case 0
        case 1
            axes(handles_axes)
            plot(t(indmax),signal(indmax),'vk','MarkerSize',10,'Linewidth',1)
            plot(t(indmin(1,:)),indmin(2,:),'^k','MarkerSize',10,'Linewidth',1)
            plot(t(gradin(1,:)),gradin(4,:),'ok','MarkerSize',10,'Linewidth',1)
    end
    
    indmaxs(j,1:length(indmax)) = indmax';
    indmins(j,1:size(indmin,2),1:2) = indmin';
    gradins(j,1:size(gradin,2),1:4) = gradin';
end

if continuous == 0
    if size(indmaxs,2) > 1
        error('Data appears to be continuous, i.e. multiple cycles.  If so, then please set the CONTINUOUS parameter to 1')
    end
end
if continuous == 1
    if size(indmaxs,2) < 2
        error('Data appears to be a single cycle.  If so, then please set the CONTINUOUS parameter to 0')
    end
end

%%%%%%%%%%%%% Determine Minima, Maxima and Max Gradients - End %%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Minima, Maxima and Max Gradients Sorting - Start %%%%%%%%%%%%%%%%%
% Ensure that similar cycles are compared
trash1 = find(gradins(1,:,1));
trash2 = find(gradins(2,:,1));

ind1 = gradins(1,trash1,1);
ind2 = gradins(2,trash2,1);

% Determine the standard deviation between the number of data points between 
% successive cycles for both signals.  Which ever has the lowest standard
% deviation will be used as the reference vector, (assuming this has the
% least noise/variabilit etc.
diff1 = ind1(1,1:end-1) - ind1(1,2:end);
diff1_std = std(diff1);

diff2 = ind2(1,1:end-1) - ind2(1,2:end);
diff2_std = std(diff2);

ind1_new = [];
ind2_new = [];

% Based on 
trash_g = [];
TTWindow = 0.25; % corresponding feet will be scanned to ensure they are within +/- 25% of the reference point

while isempty(trash_g) % In case TTWindow is too small, (e.g. if carotid/ankle waveforms; high delay)
    count = 1;
    if diff1_std < diff2_std % signal 1 = reference feet
        std_int = median(diff1); % Median number of points between reference feet
        
        l1 = length(ind1);
        l2 = length(ind2);
        for i=1:l1
            % Locate the nearest point to the reference by creating a
            % 'mask', this is done by subtracting all the non-reference indexes 
            % from the current point on the reference vector
            trash1 = abs(ones(1,l2).*ind1(i) - ind2);
            trash2 = find(trash1 > (std_int*TTWindow));
            trash1(trash2) = nan;
            if length(ind2) > length(trash2) % ensure the compared index length is greater than the mask length
                [~,I] = min(trash1);
                
                trash_m(1,count,:) = indmins(1,i,:);
                trash_m(2,count,:) = indmins(2,I,:);
                
                trash_g(1,count,:) = gradins(1,i,:);
                trash_g(2,count,:) = gradins(2,I,:);
                
                trash_n(1,count,:) = indmaxs(1,i);
                trash_n(2,count,:) = indmaxs(2,I);
                
                count = count + 1;
            end
        end
    else % signal 2 = reference feet
        std_int = median(diff2); % Median number of points between reference feet
        
        l1 = length(ind1);
        l2 = length(ind2);
        for i=1:l2
            % Locate the nearest point to the reference by creating a
            % 'mask', this is done by subtracting all the non-reference indexes 
            % from the current point on the reference vector
            trash1 = abs(ones(1,l1).*ind2(i) - ind1);
            trash2 = find(trash1 > (std_int*TTWindow));
            trash1(trash2) = nan;
            if length(ind1) > length(trash2) % ensure the compared index length is greater than the mask length
                [~,I] = min(trash1);
                
                trash_m(2,count,:) = indmins(2,i,:);
                trash_m(1,count,:) = indmins(1,I,:);
                
                trash_g(2,count,:) = gradins(2,i,:);
                trash_g(1,count,:) = gradins(1,I,:);
                
                trash_n(2,count,:) = indmaxs(2,i);
                trash_n(1,count,:) = indmaxs(1,I);
                
                count = count + 1;
            end
        end
    end
    TTWindow = TTWindow + 0.1;
end

%%%%%%%%%%%%% Minima, Maxima and Max Gradients Sorting - End %%%%%%%%%%%%%%%%%

clear gradins indmins indmaxs
gradins = trash_g;
indmins = trash_m;
indmaxs = trash_n;

switch show
    case 0
    case 1
        axes(handles_axes)
        for i=1:2
            plot(t(indmaxs(i,:)),signal1(i,indmaxs(i,:)),'vr','MarkerSize',12,'Linewidth',1)
            plot(t(indmins(i,:,1)),indmins(i,:,2),'^r','MarkerSize',12,'Linewidth',1)
            plot(t(gradins(i,:,1)),gradins(i,:,4),'or','MarkerSize',12,'Linewidth',1)
        end
end

% %Ensure that signals 1 and 2 are in the correct order : MARIE NON
% diff_ind = indmins(1,:,1) - indmins(2,:,1);
% trash = median(diff_ind);
% 
% if trash < 0
% else
%     signal1 = flipud(signal1);
%     
%     indmins = flipdim(indmins,1);
%     indmaxs = flipdim(indmaxs,1);
%     gradins = flipdim(gradins,1);
% end
signal = signal1';

function [ind_intersection] = intersection(flow)
    l = length(flow);
    m = 1;
    for k=2:l
        test = flow(k-1)*flow(k);
        if test<0
            ind_intersection(m) = k;
            m=m+1;
        end
    end

function t_pre_ejection = get_t_pre_ejection(Q)
    [c ind_max] = max(Q);
    dQ = diff(Q);
    [ind] = intersection(dQ);
    ind_zero = find(ind==ind_max);
    if ind_zero == 1
        t_pre_ejection = 1;
    else
        t_pre_ejection = ind(ind_zero-1);
    end
    
function ENd_avg = compute_ENd_avg(t)    
    c = [0.35695;-7.2266;74.249;-307.39;684.54;-856.92;571.95;-159.1]
    for k=1:length(c)
        val(k) = c(k)*t^(k-1)        
    end
    ENd_avg = sum(val);
    
function ENd_est = compute_ENd_est(EF,DBP,Pes,ENd_avg) 
    ENd_est = 0.0275 - 0.165*EF + 0.3656*DBP/Pes + 0.515*ENd_avg;

function Ees = compute_Ees(DBP,ENd_est,SBP,SV) 
    Ees = (DBP-ENd_est*SBP*0.9)/(SV*ENd_est);

function [Dis_Asc Dis_Desc Dis_Dia] = compute_distensibility(PP,areamm)
    %- Ascending
    Amin_asc = min(areamm(:,1)); Amax_asc = max(areamm(:,1));
    Dis_Asc = (Amax_asc-Amin_asc)/Amin_asc/PP;
    %- Descending
    Amin_desc = min(areamm(:,2)); Amax_desc = max(areamm(:,2));
    Dis_Desc = (Amax_desc-Amin_desc)/Amin_desc/PP;
    %- Diaphram
    Amin_dia = min(areamm(:,3)); Amax_dia = max(areamm(:,3));
    Dis_Dia = (Amax_dia-Amin_dia)/Amin_dia/PP;

function Generate_spreadsheet_button_Callback(hObject, eventdata, handles)
global transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia length_asc_desc length_desc_dia length_asc_dia ...
        transit_time path_length PWV Tasc Tdesc Tdia Qasc Qdesc Qdia Aasc Adesc Adia PWV_asc_desc PWV_desc_dia PWV_asc_dia ...
        info Dis_Asc Dis_Desc Dis_Dia Tasc Tdesc Tdia Qasc Qdesc Qdia Aasc Adesc Adia EDV EF PP DBP SBP Pes ENd_est Ees tNd ENd_avg
EF1 = getappdata(0,'EF1'); HR = getappdata(0,'HR'); SV = getappdata(0,'SV'); CO = getappdata(0,'CO');
Dis_Asc = getappdata(0,'Dis_Asc'); Dis_Desc = getappdata(0,'Dis_Desc'); Dis_Dia = getappdata(0,'Dis_Dia');
info = getappdata(0,'info');
path_flow = getappdata(0,'path_flow')
path_origin = cd;
if isempty(Tdia)
    Time_s = Tasc; Flow_Ascending_ml_s = Qasc; Flow_Descending_ml_s = Qdesc;
    Area_Ascending_mm2 = Aasc; Area_Descending_mm2 = Adesc;
    Length_Ascending_Descending_m = length_asc_desc; TT_Ascending_Descending_s = transit_time_asc_desc; PWV_Ascending_Descending_m_s = PWV_asc_desc;
    EF1_pct = EF1; SV_ml = SV; HR_bpm = HR; CO_m3_s = CO;
    cd(path_flow);
    [file,path,indx] = uiputfile('');
    generate_spreadsheet_1_parts(file,path,info,EF1,SV,HR,CO,transit_time_asc_desc,Dis_Asc,Dis_Desc,Dis_Dia,Tasc,Qasc,Qdesc,Qdia,Aasc,Adesc,Adia,EDV,EF,PP,DBP,SBP,Pes,ENd_est,Ees,tNd,ENd_avg)
else    
    Time_s = Tasc; Flow_Ascending_ml_s = Qasc; Flow_Descending_ml_s = Qdesc; Flow_Diaphragm_ml_s = Qdia;
    Area_Ascending_mm2 = Aasc; Area_Descending_mm2 = Adesc; Area_Diaphragm_mm2 = Adia; 
    TT_Ascending_Descending_s = transit_time_asc_desc; TT_Descending_Diaphragm_s = transit_time_desc_dia; TT_Ascending_Diaphragm_s = transit_time_asc_dia;
    PWV_Ascending_Descending_m_s = PWV_asc_desc; PWV_Descending_Diaphragm_m_s = PWV_desc_dia; PWV_Ascending_Diaphragm_m_s = PWV_asc_dia;
    EF1_pct = EF1; SV_ml = SV; HR_bpm = HR; CO_m3_s = CO;
    cd(path_flow);
    [file,path,indx] = uiputfile('.xls');
    generate_spreadsheet_3_parts(file,path,info,EF1,SV,HR,CO,transit_time_asc_desc,transit_time_desc_dia,transit_time_asc_dia,Dis_Asc,Dis_Desc,Dis_Dia,Tasc,Qasc,Qdesc,Qdia,Aasc,Adesc,Adia,EDV,EF,PP,DBP,SBP,Pes,ENd_est,Ees,tNd,ENd_avg)
end
cd(path_origin)

function [] = generate_spreadsheet_1_part(filename,path,info,EF1,SV,HR,CO,transit_time_asc_desc,Dis_Asc,Dis_Desc,Tasc,Qasc,Qdesc,Aasc,Adesc,EDV,EF,PP,DBP,SBP,Pes,ENd_est,Ees,tNd,ENd_avg)

cd(path)

fill = {'','','','','',''};

%% PATIENT INFO
%- Check fields exist
if isfield(info,'PatientName')
    if isfield(info,'PatientName.GivenName')
        Field_Name = [info.PatientName.GivenName,' ',info.PatientName.FamilyName];
    else
        Field_Name = [info.PatientName.FamilyName];
    end
else
    Field_Name = '';
end
if isfield(info,'PatientID')
    Field_PatientID = [info.PatientID];
else
    Field_PatientID = '';
end
if isfield(info,'FileModDate')
    Field_ScanDate = [info.FileModDate];
else
    Field_ScanDate = '';
end
if isfield(info,'PatientBirthDate')
    if isempty(info.PatientBirthDate)==0
        Field_DOB = [info.PatientBirthDate(end-1:end),'-',info.PatientBirthDate(end-3:end-2),'-',info.PatientBirthDate(1:4)];
    else 
        Field_DOB = '';
    end
else
    Field_DOB = '';
end
if isfield(info,'AccessionNumber')
    Field_AccessionNumber = [info.AccessionNumber];
else
    Field_AccessionNumber = '';
end

line_Patient_Info = [{'%%% PATIENT INFO','','','','',''};
       {'Name:',Field_Name,'','','',''};
       {'Patient ID:',Field_PatientID,'','','',''};
       {'Scan date:',Field_ScanDate,'','','',''};
       {'DOB:',Field_DOB,'','','',''};
       {'Accession number:',Field_AccessionNumber,'','','',''};
       ];

%% Transit time
line_Transit_Time = [{'%%% TRANSIT TIME','','','','',''};
       {'Transit time ascending-descending aorta:',transit_time_asc_desc,'s','','',''};
       {'Transit time descending-diaphragm aorta:',transit_time_desc_dia,'s','','',''};
       {'Transit time ascending-diaphragm aorta:',transit_time_asc_dia,'s','','',''};
       ];  
   
   %% FLOW ANALYSIS
line_Flow_Analysis = [{'%%% FLOW ANALYSIS','','','','',''};
       {'Stroke volume:',round(100*SV)/100,'ml','','',''};
       {'Heart rate:',round(100*HR)/100,'bpm','','',''};
       {'Cardiac output:',round(100*CO)/100,'L/min','','',''};
       {'EF1:',round(100*EF1)/100,'%','','',''};
       ];
   
%% INPUTS
line_Inputs = [{'%%% INPUTS','','','','',''};
       {'DBP:',DBP,'mmHg','','',''};
       {'SBP:',SBP,'mmHg','','',''};
       {'PP:',PP,'mmHg','','',''};
       {'Pes:',Pes,'mmHg','','',''};
       {'EDV:',round(EDV),'ml','','',''};
       {'EF:',round(EF*1e4),'%','','',''};
       ];   

%% LV Elastance
line_Elastance = [{'%%% ELASTANCE','','','','',''};
       {'Elastance at onset of ejection (ENd):',round(ENd_est*100)/100,'mmHg/ml','','',''};
       {'Elastance at end systole (Ees):',round(Ees*100)/100,'mmHg/ml','','',''};
       ];   
   
%% Arterial Distensibility
line_Distensibility = [{'%%% ARTERIAL DISTENSIBILITY','','','','',''};
       {'Ascending aorta:',Dis_Asc,'/mmHg','','',''};
       {'Descending aorta:',Dis_Desc,'/mmHg','','',''};
       {'Diaphragm:',Dis_Dia,'/mmHg','','',''};
       ];   
   
%% Flow waves
filler = strings(length(Qasc),1)
line_Flow_Waves = [["%%% FLOW WAVES","","","","",""];["Time (s)","Ascending aorta (ml/s)","Descending aorta (ml/s)","","",""];...
    [Tasc/1e3,Qasc,Qdesc,filler,filler,filler]]

%% Luminal variations waves
filler = strings(length(Aasc),1)
line_Luminal_Variations = [["%%% Luminal variations","","","","",""];["Time (s)","Ascending aorta (mm2)","Descending aorta (mm2)","","",""];...
    [Tasc/1e3,Aasc,Adesc,filler,filler,filler]]

%% ASSEMBLE TABLE
T = [line_Patient_Info ; fill ; line_Transit_Time ; fill ; line_Flow_Analysis ; fill ; line_Inputs ; fill ; line_Elastance ; fill ; line_Distensibility ; fill ; cellstr(line_Flow_Waves) ; fill ; cellstr(line_Luminal_Variations)];

%% WRITE FILE
writematrix(T,[path filename]);


function [] = generate_spreadsheet_3_parts(filename,path,info,EF1,SV,HR,CO,transit_time_asc_desc,transit_time_desc_dia,transit_time_asc_dia,Dis_Asc,Dis_Desc,Dis_Dia,Tasc,Qasc,Qdesc,Qdia,Aasc,Adesc,Adia,EDV,EF,PP,DBP,SBP,Pes,ENd_est,Ees,tNd,ENd_avg)

cd(path)

fill = {'','','','','',''};

%% PATIENT INFO
%- Check fields exist
if isfield(info,'PatientName')
    if isfield(info,'PatientName.GivenName')
        Field_Name = [info.PatientName.GivenName,' ',info.PatientName.FamilyName];
    else
        Field_Name = [info.PatientName.FamilyName];
    end
else
    Field_Name = '';
end
if isfield(info,'PatientID')
    Field_PatientID = [info.PatientID];
else
    Field_PatientID = '';
end
if isfield(info,'FileModDate')
    Field_ScanDate = [info.FileModDate];
else
    Field_ScanDate = '';
end
if isfield(info,'PatientBirthDate')
    if isempty(info.PatientBirthDate)==0
        Field_DOB = [info.PatientBirthDate(end-1:end),'-',info.PatientBirthDate(end-3:end-2),'-',info.PatientBirthDate(1:4)];
    else 
        Field_DOB = '';
    end
else
    Field_DOB = '';
end
if isfield(info,'AccessionNumber')
    Field_AccessionNumber = [info.AccessionNumber];
else
    Field_AccessionNumber = '';
end

line_Patient_Info = [{'%%% PATIENT INFO','','','','',''};
       {'Name:',Field_Name,'','','',''};
       {'Patient ID:',Field_PatientID,'','','',''};
       {'Scan date:',Field_ScanDate,'','','',''};
       {'DOB:',Field_DOB,'','','',''};
       {'Accession number:',Field_AccessionNumber,'','','',''};
       ];

%% Transit time
line_Transit_Time = [{'%%% TRANSIT TIME','','','','',''};
       {'Transit time ascending-descending aorta:',transit_time_asc_desc,'s','','',''};
       {'Transit time descending-diaphragm aorta:',transit_time_desc_dia,'s','','',''};
       {'Transit time ascending-diaphragm aorta:',transit_time_asc_dia,'s','','',''};
       ];  
   
   %% FLOW ANALYSIS
line_Flow_Analysis = [{'%%% FLOW ANALYSIS','','','','',''};
       {'Stroke volume:',round(SV),'ml','','',''};
       {'Heart rate:',round(HR),'bpm','','',''};
       {'Cardiac output:',round(100*CO)/100,'L/min','','',''};
       {'EF1:',round(EF1),'%','','',''};
       ];
   
%% INPUTS
line_Inputs = [{'%%% INPUTS','','','','',''};
       {'DBP:',round(DBP),'mmHg','','',''};
       {'SBP:',round(SBP),'mmHg','','',''};
       {'PP:',round(PP),'mmHg','','',''};
       {'Pes:',round(Pes),'mmHg','','',''};
       {'EDV:',round(EDV),'ml','','',''};
       {'EF:',round(EF),'%','','',''};
       ];   

%% LV Elastance
line_Elastance = [{'%%% ELASTANCE','','','','',''};
       {'Elastance at onset of ejection (ENd):',round(ENd_est*100)/100,'mmHg/ml','','',''};
       {'Elastance at end systole (Ees):',round(Ees*100)/100,'mmHg/ml','','',''};
       ];   
   
%% Arterial Distensibility
line_Distensibility = [{'%%% ARTERIAL DISTENSIBILITY','','','','',''};
       {'Ascending aorta:',Dis_Asc,'/mmHg','','',''};
       {'Descending aorta:',Dis_Desc,'/mmHg','','',''};
       {'Diaphragm:',Dis_Dia,'/mmHg','','',''};
       ];   
   
%% Flow waves
filler = strings(length(Qasc),1)
line_Flow_Waves = [["%%% FLOW WAVES","","","","",""];["Time (s)","Ascending aorta (ml/s)","Descending aorta (ml/s)","Diaphragm (ml/s)","",""];...
    [Tasc/1e3,Qasc,Qdesc,Qdia,filler,filler]]

%% Luminal variations waves
filler = strings(length(Aasc),1)
line_Luminal_Variations = [["%%% Luminal variations","","","","",""];["Time (s)","Ascending aorta (mm2)","Descending aorta (mm2)","Diaphragm (mm2)","",""];...
    [Tasc/1e3,Aasc,Adesc,Adia,filler,filler]]

%% ASSEMBLE TABLE
T = [line_Patient_Info ; fill ; line_Transit_Time ; fill ; line_Flow_Analysis ; fill ; line_Inputs ; fill ; line_Elastance ; fill ; line_Distensibility ; fill ; cellstr(line_Flow_Waves) ; fill ; cellstr(line_Luminal_Variations)];

%% WRITE FILE
writecell(T,[path filename]);
