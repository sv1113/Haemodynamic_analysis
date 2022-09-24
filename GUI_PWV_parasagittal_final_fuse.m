function varargout = GUI_PWV_parasagittal_final(varargin)
% GUI_PWV_PARASAGITTAL_FINAL MATLAB code for GUI_PWV_parasagittal_final.fig
%      GUI_PWV_PARASAGITTAL_FINAL, by itself, creates a new GUI_PWV_PARASAGITTAL_FINAL or raises the existing
%      singleton*.
%
%      H = GUI_PWV_PARASAGITTAL_FINAL returns the handle to a new GUI_PWV_PARASAGITTAL_FINAL or the handle to
%      the existing singleton*.
%
%      GUI_PWV_PARASAGITTAL_FINAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PWV_PARASAGITTAL_FINAL.M with the given input arguments.
%
%      GUI_PWV_PARASAGITTAL_FINAL('Property','Value',...) creates a new GUI_PWV_PARASAGITTAL_FINAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_PWV_parasagittal_final_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_PWV_parasagittal_final_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_PWV_parasagittal_final

% Last Modified by GUIDE v2.5 09-Mar-2020 15:52:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_PWV_parasagittal_final_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_PWV_parasagittal_final_OutputFcn, ...
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

% --- Executes just before GUI_PWV_parasagittal_final is made visible.
function GUI_PWV_parasagittal_final_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_PWV_parasagittal_final (see VARARGIN)

% Choose default command line output for GUI_PWV_parasagittal_final
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_PWV_parasagittal_final wait for user response (see UIRESUME)
global transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia length_asc_desc length_desc_dia length_asc_dia PWV_asc_desc PWV_desc_dia PWV_asc_dia
transit_time_asc_desc = []; transit_time_desc_dia = []; transit_time_asc_dia = []; 
length_asc_desc = []; length_desc_dia = []; length_asc_dia = []; 
PWV_asc_desc = [];  PWV_desc_dia = []; PWV_asc_dia = []; 
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_PWV_parasagittal_final_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global x_centreline y_centreline PixelSpacing transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia length_asc_desc length_desc_dia length_asc_dia PWV_asc_desc PWV_desc_dia PWV_asc_dia
setappdata(0,'x_centreline',x_centreline)
setappdata(0,'y_centreline',y_centreline)
setappdata(0,'PixelSpacing',PixelSpacing)
setappdata(0,'transit_time_asc_desc',transit_time_asc_desc)
setappdata(0,'transit_time_desc_dia',transit_time_desc_dia)
setappdata(0,'transit_time_asc_dia',transit_time_asc_dia)
% Get default command line output from handles structure
handles.output = [{x_centreline} {y_centreline} PixelSpacing]
varargout{1} = handles.output;
if (~isempty(x_centreline)) && (~isempty(y_centreline)) && (~isempty(PixelSpacing))
    uiresume
%     close
end
%- Store variable in handle
guidata(hObject, handles);

% --- Executes on button press in Load_scans_button.
function Load_scans_button_Callback(hObject, eventdata, handles)
global info Xf dname_folder axsize figsize x_up y_up x_down y_down transit_time path_length PWV current_folder
x_up = []; y_up = []; x_down = []; y_down = []; transit_time = []; path_length = []; PWV = [];
current_folder = cd
dname = uigetdir(current_folder,'Select folder containing parasagittal scans');
cd(dname)
files = dir('*.dcm');
[Xf info] = fuse_images(files);
axes(handles.axes1)
imshow(Xf);
axsize=get(gca,'position');
figsize=get(gcf,'position'); 
cd(current_folder)
axes(handles.axes2)
guidata(hObject, handles);

% --- Executes on button press in Reference_slice_transit_time_button.
function Reference_slice_transit_time_button_Callback(hObject, eventdata, handles)
global z_ref Xf info dname_folder z_ref_matrix current_folder PixelSpacing z_ref_diaphragm_matrix transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia length_asc_desc length_desc_dia length_asc_dia
current_folder = cd;
[filename1,filepath1]=uigetfile({'*.*','All Files'},...
  'Select Flow File');
global transit_time path_length PWV

%- Display reference line
cd(filepath1)
files = dir('*.dcm');
info_axial = dicominfo(files(1).name);
% cd(dname_folder)
z_ref = info_axial.SliceLocation;
PatientPosition = info.ImagePositionPatient;
PixelSpacing = info.PixelSpacing(2);
z_loc = PatientPosition(3) - [0:size(Xf,1)-1].*PixelSpacing;
[c z_ref_matrix] = min(abs(z_ref - z_loc));
Xf(z_ref_matrix,:,1)=0; Xf(z_ref_matrix,:,2)=0; Xf(z_ref_matrix,:,3)=1;
%- Artificially set diaphrahgm limit
nb_shift = 70;
z_ref_diaphragm_matrix = z_ref_matrix+nb_shift;
Xf(z_ref_diaphragm_matrix,:,1)=0; Xf(z_ref_diaphragm_matrix+nb_shift,:,2)=0; Xf(z_ref_diaphragm_matrix+nb_shift,:,3)=1;

%- Find transit time
[Tasc Tdesc Tdia Qasc Qdesc Qdia] = extract_Q_wave([filepath1 filename1])
cd(current_folder)
if Tasc(1) == 0
    freq = 1e3/Tasc(2);
else
    freq = 1e3/Tasc(1);
end
%- Open GUI for calculating transit time
if isempty(Tdia)
    
else
    flow = [Qasc Qdesc Qdia];
    setappdata(0,'flow',flow)
    setappdata(0,'freq',freq)
    GUI_Transit_Time_3_signals
    %- Display TT ascending descending
    transit_time_asc_desc = getappdata(0,'transit_time_asc_desc');
    str_TT_asc_desc = sprintf('%.3f',transit_time_asc_desc)
    set(handles.str_TT_asc_desc,'String',str_TT_asc_desc);
    if ~isempty(length_asc_desc)==1
        PWV_asc_desc = length_asc_desc*1e-3/transit_time_asc_desc;
        PWV_asc_desc_str = sprintf('%.1f', PWV_asc_desc);
        set(handles.str_PWV_asc_desc, 'String', PWV_asc_desc_str);
    end
    %- Display TT descending diaphragm
    transit_time_desc_dia = getappdata(0,'transit_time_desc_dia');
    str_TT_desc_dia = sprintf('%.3f',transit_time_desc_dia)
    set(handles.str_TT_desc_dia,'String',str_TT_desc_dia);
    if ~isempty(length_desc_dia)==1
        PWV_desc_dia = length_desc_dia*1e-3/transit_time_desc_dia;
        PWV_desc_dia_str = sprintf('%.1f', PWV_desc_dia);
        set(handles.str_PWV_desc_dia, 'String', PWV_desc_dia_str);
    end
    %- Display TT ascending diaphragm
    transit_time_asc_dia = getappdata(0,'transit_time_asc_dia');
    str_TT_asc_dia = sprintf('%.3f',transit_time_asc_dia)
    set(handles.str_TT_asc_dia,'String',str_TT_asc_dia);
    if ~isempty(length_asc_dia)==1
        PWV_asc_dia = length_asc_dia*1e-3/transit_time_asc_dia;
        PWV_asc_dia_str = sprintf('%.1f', PWV_asc_dia);
        set(handles.str_PWV_asc_desc, 'String', PWV_asc_dia_str);
    end
end

axes(handles.axes1)
imshow(Xf);
cd(current_folder)

%- Store variable in handle
guidata(hObject, handles);



% --- Executes on button press in Reset_button.
function Reset_button_Callback(hObject, eventdata, handles)
global transit_time path_length PWV
axes(handles.axes1); cla reset;
axes(handles.axes2); cla reset;
set(handles.str_length_asc_desc,'String','');set(handles.str_length_desc_dia,'String','');set(handles.str_length_asc_dia,'String','');
set(handles.str_PWV_asc_desc,'String','');set(handles.str_PWV_desc_dia,'String','');set(handles.str_PWV_asc_dia,'String','');
set(handles.str_TT_asc_desc,'String','');set(handles.str_TT_desc_dia,'String','');set(handles.str_TT_asc_dia,'String','');
clear transit_time path_length PWV
%- Store variable in handle
guidata(hObject, handles);


% --- Executes on slider movement.
function Threshold_asc_slider_Callback(hObject, eventdata, handles)
global percentage_threshold Xf Aorta Xgray z_ref z_ref_matrix x_up y_up x_down y_down x_centre x1 X_segmented_aorta
val = get(hObject,'Value');
val_min = get(hObject,'Min'); val_max = get(hObject,'Max');
percentage_threshold = (val - val_min)/(val_max - val_min);
Xf = getimage(handles.axes1);
Xgray = rgb2gray(Xf);
threshold = max(Xgray(:))*percentage_threshold;
[Aorta xline yline] = extract_aorta(Xgray,threshold);
X_segmented_aorta = imfuse(Xf,Aorta,'blend');
axes(handles.axes2)
imshow(Xf)
[x_up y_up x_down y_down x_centre x1] = coordinate_boundaries(xline, yline)
hold all, plot(x_up,y_up,'b'), plot(x_down,y_down,'b')
%- Store variable in handle
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Threshold_asc_slider_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in Upper_boundary_asc_button.
function Upper_boundary_asc_button_Callback(hObject, eventdata, handles)
global Xf x_up y_up x_down y_down
axes(handles.axes2)
imshow(Xf);
if isempty(x_down)==0 && isempty(y_down)==0
    plot(x_down,y_down,'b')
end
% datacursormode on
set(gcf,'windowbuttondownfcn',{@starttrack});
waitfor(0,'userdata')
upper_boundary=get(0,'userdata');
x_up = upper_boundary(:,1); 
y_up = handles.axes2.YLim(2)-upper_boundary(:,2); 
hold on,plot(x_up,y_up,'b'),ylim(handles.axes2.YLim),xlim(handles.axes2.XLim),set(gca,'Position',handles.axes2.Position)
%- Store variable in handle
guidata(hObject, handles);

% --- Executes on button press in Lower_boundary_asc_button.
function Lower_boundary_asc_button_Callback(hObject, eventdata, handles)
global Xf x_up y_up x_down y_down
axes(handles.axes2)
imshow(Xf);
if isempty(x_up)==0 && isempty(y_up)==0
    plot(x_up,y_up,'b')
end
% datacursormode on
set(gcf,'windowbuttondownfcn',{@starttrack});
waitfor(0,'userdata')
upper_boundary=get(0,'userdata');
x_down = upper_boundary(:,1); 
y_down = handles.axes2.YLim(2)-upper_boundary(:,2); 
hold on,plot(x_down,y_down,'b'),ylim(handles.axes2.YLim),xlim(handles.axes2.XLim),set(gca,'Position',handles.axes2.Position)
%- Store variable in handle
guidata(hObject, handles);

% --- Executes on button press in Algorithm_asc_button.
function Algorithm_asc_button_Callback(hObject, eventdata, handles)
global x_up y_up x_down y_down x_centre x1 z_ref_matrix Xf info aortic_length x_centreline y_centreline transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia length_asc_desc length_desc_dia length_asc_dia PWV_asc_desc PWV_desc_dia PWV_asc_dia
x_centre = (y_down(1)+y_down(end))/2;
flag_plot = 1;
[centreline] = find_centreline(x_down,y_down,x_up,y_up,x_centre,z_ref_matrix,flag_plot)
x_centreline = centreline(:,1);
y_centreline = -centreline(:,2) + z_ref_matrix;
%- Show centreline
axes(handles.axes2), hold all
imshow(Xf),plot(x_up,y_up,'b'),plot(x_down,y_down,'b'),plot(x_centreline,y_centreline,'r','linewidth',2)
length_asc_desc = calculate_length(centreline(:,1),centreline(:,2),info)
str_length_asc_desc = sprintf('%.3f', length_asc_desc/1e3);
set(handles.str_length_asc_desc, 'String', str_length_asc_desc);
if ~isempty(transit_time_asc_desc)==1
    PWV_asc_desc = length_asc_desc*1e-3/transit_time_asc_desc;
    pwv_asc_desc_str = sprintf('%.1f', PWV_asc_desc);
    set(handles.str_PWV_asc_desc, 'String', pwv_asc_desc_str);
    uiresume
end
%- Store variable in handle
guidata(hObject, handles);

% --- Executes on button press in Manual_tracing_asc_button.
function Manual_tracing_asc_button_Callback(hObject, eventdata, handles)
global Xf x_up y_up x_down y_down   info z_ref_matrix transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia length_asc_desc length_desc_dia length_asc_dia PWV_asc_desc PWV_desc_dia PWV_asc_dia
axes(handles.axes2)
imshow(Xf);
if isempty(x_down)==0 && isempty(y_down)==0
    plot(x_down,y_down,'b')
end
if isempty(x_up)==0 && isempty(y_up)==0
    plot(x_up,y_up,'b')
end
% datacursormode on
set(gcf,'windowbuttondownfcn',{@starttrack});
waitfor(0,'userdata')
centreline=get(0,'userdata');
x_centreline_manual = centreline(:,1);
y_centreline_manual = handles.axes2.YLim(2)-centreline(:,2); 

%- Delete what is below reference line (z_ref_matrix)
ind = find(y_centreline_manual>z_ref_matrix);
x_centreline_manual(ind) = [];
y_centreline_manual(ind) = [];
x_centreline_manual = [x_centreline_manual(1);x_centreline_manual;x_centreline_manual(end)];
y_centreline_manual = [z_ref_matrix;y_centreline_manual;z_ref_matrix];

hold on,plot(x_centreline_manual,y_centreline_manual,'r','linewidth',2),
ylim(handles.axes2.YLim),xlim(handles.axes2.XLim),
set(gca,'Position',handles.axes2.Position)

setappdata(0,'x_centreline',x_centreline)
setappdata(0,'y_centreline',y_centreline)
setappdata(0,'PixelSpacing',PixelSpacing)

length_asc_desc = calculate_length(x_centreline_manual,y_centreline_manual,info)
str_length_asc_desc = sprintf('%.3f', length_asc_desc/1e3);
set(handles.str_length_asc_desc, 'String', str_length_asc_desc);
if ~isempty(transit_time_asc_desc)==1
    PWV_asc_desc = length_asc_desc*1e-3/transit_time_asc_desc;
    pwv_asc_desc_str = sprintf('%.1f', PWV_asc_desc);
    set(handles.str_PWV_asc_desc, 'String', pwv_asc_desc_str);
    uiresume
end
%- Store variable in handle
guidata(hObject, handles);


function str_TT_asc_desc_Callback(hObject, eventdata, handles)
global transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia length_asc_desc length_desc_dia length_asc_dia PWV_asc_desc PWV_desc_dia PWV_asc_dia
transit_time_asc_desc = str2double(get(hObject,'String'));
if ~isempty(length_asc_desc)==1
    PWV_asc_desc = length_asc_desc*1e-3/transit_time_asc_desc;
    pwv_asc_desc_str = sprintf('%.1f', PWV_asc_desc);
    set(handles.str_PWV_asc_desc, 'String', pwv_asc_desc_str);
    uiresume
end
%- Store variable in handle
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function str_TT_asc_desc_CreateFcn(hObject, eventdata, handles)

% --- Executes on slider movement to segment the descending - diaphragm part.
function Threshold_diaphragm_slider_Callback(hObject, eventdata, handles)
global percentage_threshold Xf Aorta Xgray z_ref z_ref_matrix X_segmented_aorta z_ref_diaphragm_matrix x_left y_left x_right y_right
val = get(hObject,'Value');
val_min = get(hObject,'Min'); val_max = get(hObject,'Max');
percentage_threshold = (val - val_min)/(val_max - val_min);
Xf = getimage(handles.axes1);
Xgray = rgb2gray(Xf);
threshold = max(Xgray(:))*percentage_threshold;
[Lower_Aorta xline_lower yline_lower] = extract_diaphragm_aorta(Xgray,threshold);
Aorta = [zeros(z_ref_matrix-1,size(Aorta,2)) ; Lower_Aorta ; zeros(size(Xf,1)-z_ref_diaphragm_matrix,size(Aorta,2))]
X_segmented_aorta = imfuse(Xf,Aorta,'blend');
axes(handles.axes2)
imshow(Xf)
[x_left y_left x_right y_right] = coordinate_boundaries_lower(xline_lower, yline_lower,z_ref_matrix,z_ref_diaphragm_matrix)
hold all, plot(x_left,y_left,'b'), plot(x_right,y_right,'b')
%- Store variable in handle
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Threshold_diaphragm_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in Left_boundary_diaphragm_button.
function Left_boundary_diaphragm_button_Callback(hObject, eventdata, handles)
global Xf x_left y_left x_right y_right info z_ref_matrix transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia length_asc_desc length_desc_dia length_atransit_time_asc_desc transit_time_desc_dia transit_time_asc_dia length_asc_desc length_desc_dia length_asc_dia PWV_asc_desc PWV_desc_dia PWV_asc_diasc_dia PWV_asc_desc PWV_desc_dia PWV_asc_dia
axes(handles.axes2)
imshow(Xf); hold all;
if isempty(x_right)==0 && isempty(y_right)==0
    plot(x_right,y_right,'b')
end
% datacursormode on
set(gcf,'windowbuttondownfcn',{@starttrack});
waitfor(0,'userdata')
upper_boundary=get(0,'userdata');
x_left = upper_boundary(:,1); 
y_left = handles.axes2.YLim(2)-upper_boundary(:,2); 
hold on,plot(x_left,y_left,'b'),ylim(handles.axes2.YLim),xlim(handles.axes2.XLim),set(gca,'Position',handles.axes2.Position)
%- Store variable in handle
guidata(hObject, handles);

% --- Executes on button press in Right_boundary_diaphragm_button.
function Right_boundary_diaphragm_button_Callback(hObject, eventdata, handles)
global Xf x_left y_left x_right y_right info z_ref_matrix transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia length_asc_desc length_desc_dia length_asc_dia PWV_asc_desc PWV_desc_dia PWV_asc_dia
axes(handles.axes2)
imshow(Xf); hold all;
if isempty(x_left)==0 && isempty(y_left)==0
    plot(x_left,y_left,'b')
end
% datacursormode on
set(gcf,'windowbuttondownfcn',{@starttrack});
waitfor(0,'userdata')
upper_boundary=get(0,'userdata');
x_right = upper_boundary(:,1); 
y_right = handles.axes2.YLim(2)-upper_boundary(:,2); 
hold on,plot(x_right,y_right,'b'),ylim(handles.axes2.YLim),xlim(handles.axes2.XLim),set(gca,'Position',handles.axes2.Position)
%- Store variable in handle
guidata(hObject, handles);

% --- Executes on button press in Algorithm_diaphragm_button.
function Algorithm_diaphragm_button_Callback(hObject, eventdata, handles)
global x_left y_left x_right y_right Xf z_ref_matrix z_ref_diaphragm_matrix info transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia length_asc_desc length_desc_dia length_asc_dia PWV_asc_desc PWV_desc_dia PWV_asc_dia
flag_plot = 1;
[centreline] = find_centreline_diaphragm(x_left,y_left,x_right,y_right,flag_plot);
axes(handles.axes2)
imshow(Xf), hold all, plot(x_left,y_left,'b'), plot(x_right,y_right,'b'), plot(centreline(1,:),centreline(2,:),'r')
length_desc_dia = calculate_length(centreline(:,1),centreline(:,2),info)
text_length_desc_dia = sprintf('%.3f', length_desc_dia/1e3);
set(handles.str_length_desc_dia, 'String', text_length_desc_dia);
if ~isempty(transit_time_desc_dia)==1
    PWV_desc_dia = length_desc_dia*1e-3/transit_time_desc_dia;
    pwv_desc_dia_str = sprintf('%.1f', PWV_desc_dia);
    set(handles.str_PWV_desc_dia, 'String', pwv_desc_dia_str);
    if ~isempty(PWV_asc_desc)==1
        length_asc_dia = length_asc_desc + length_desc_dia;
        length_asc_dia_str = sprintf('%.1f', length_asc_dia*1e-3);
        set(handles.str_length_asc_dia, 'String', length_asc_dia_str);
        PWV_asc_dia = length_asc_dia*1e-3/transit_time_asc_dia;
        pwv_asc_dia_str = sprintf('%.1f', PWV_asc_dia);
        set(handles.str_PWV_asc_dia, 'String', pwv_asc_dia_str);
    end
    uiresume
end
%- Store variable in handle
guidata(hObject, handles);

% --- Executes on button press in Manual_tracing_diaphragm_button.
function Manual_tracing_diaphragm_button_Callback(hObject, eventdata, handles)

function str_TT_desc_dia_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function str_TT_desc_dia_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% ADDITIONAL FUNCTIONS
function [Xf info] = fuse_images(files)
    for k=1:round(size(files,1))
        X = dicomread(files(k).name);
        info = dicominfo(files(k).name);
        if k==1
            Xf = zeros(size(X,1),size(X,2));
        end 
        X = double(X);
    %     Xgray = rgb2gray(X);
    
        cmap = copper(256);
        recolored_image = ind2rgb(X, cmap);
        Xgray = rgb2gray(recolored_image);
        BW = edge(Xgray,'Canny');
        Xf = imfuse(Xf,recolored_image,'blend');
    %     subplot(6,5,k),imshow(Xf);
    end

function  [Aorta xline yline] = extract_aorta(Xgray,threshold)
    for i=1:size(Xgray,1)
        for j=1:size(Xgray,2)
            if Xgray(i,j)<threshold
                Xgray(i,j)=0;
            end
        end
    end
%     figure(),imshow(Xgray)            
    %% Find the aorta (largest region)
    global Aorta
    [B,L] = bwboundaries(Xgray,'holes');
    %     imshow(label2rgb(L, @jet, [.5 .5 .5]))
    for k=1:size(B,1)
        l(k) = size(B{k},1);    
    end
    [c ind] = max(l);
    temp = zeros(size(Xgray,1),size(Xgray,2));  
    for k=1:size(B{ind},1)
        temp(B{ind}(k,1),B{ind}(k,2))=1;
    end
    Aorta = temp;
    xline = B{ind}(:,1);
    yline = B{ind}(:,2);

function [Lower_Aorta xline yline] = extract_diaphragm_aorta(Xgray,threshold);
    for i=1:size(Xgray,1)
        for j=1:size(Xgray,2)
            if Xgray(i,j)<threshold
                Xgray(i,j)=0;
            end
        end
    end
    %% Find the aorta (largest region)
    global Lower_Aorta
    [B,L] = bwboundaries(Xgray,'holes');
    %     imshow(label2rgb(L, @jet, [.5 .5 .5]))
    for k=1:size(B,1)
        ymin(k) = min(B{k}(:,1));    
        ymax(k) = max(B{k}(:,1));    
        spread(k) = ymax(k) - ymin(k);
    end
    [c ind] = max(spread);
    temp = zeros(size(Xgray,1),size(Xgray,2));  
    for k=1:size(B{ind},1)
        temp(B{ind}(k,1),B{ind}(k,2))=1;
    end
    Lower_Aorta = temp;
    xline = B{ind}(:,1);
    yline = B{ind}(:,2);
    dx = diff(xline);
    ind = find(dx==0);
    xline(ind) = [];
    yline(ind) = [];
    
function [x_up y_up x_down y_down x_centre x1] = coordinate_boundaries(xline, yline)
%- Find baselines
x = yline; y=xline;
xmin = max(y);
indx_base = find(y==xmin);
indx_diff = diff(indx_base);
[pks,locs] = findpeaks(indx_diff);
locs = [1 locs' length(indx_base)];
m=1;
for k=2:length(locs)
    indBaseline{m} = indx_base(locs(k-1)+1:locs(k));
    m=m+1;
end
indBaseline{1} = [indx_base(1);indBaseline{1}];
l = cellfun(@numel, indBaseline);
[b ind] = sort(l);
Baseline1 = indBaseline{ind(end-1)};
Baseline2 = indBaseline{ind(end)};
outer_ind = [Baseline1(end) Baseline2(1) Baseline1(1) Baseline2(end)];
i1 = 1; i2 = 2;
x_down = x(outer_ind(i1):outer_ind(i2))
y_down = y(outer_ind(i1):outer_ind(i2))
% close all,figure(),hold all,plot(x_down,y_down),
% plot(xline(outer_ind(i1)),yline(outer_ind(i1)),'*r')
i1 = 3; i2 = 4;
x_up = x(1:outer_ind(i1))
y_up = y(1:outer_ind(i1))
% close all,figure(),hold all,plot(x_up,y_up), plot(x_down,y_down)
% plot(xline(outer_ind(i1)),yline(outer_ind(i1)),'*r')

%- X1
Base = [x(Baseline1);x(Baseline2)];
x1 = min(Base);

%- First boundary
x_centre = (y_down(1)+y_down(end))/2;

%- Dfunction dataout = scaledata(datain,minval,maxval)
    %- delete connection to the carotid 
    [c imin] = min(y_up);
    l = round(length(x_up)/8);
    ind_delete = [imin-l:imin+l];
    y_ok1_temp = y_up;x_ok1_temp = x_up;
    x_ok1_temp(ind_delete)=[];
    y_ok1_temp(ind_delete)=[];
    
    %- Interpolate 2nd order polynomial
    [cmin imin] = min(y_ok1_temp)
    x_pts = [x_ok1_temp(1) x_ok1_temp(imin) x_ok1_temp(end)]
    y_pts = [y_ok1_temp(1) y_ok1_temp(imin) y_ok1_temp(end)]
    xq = [x_ok1_temp(1):x_ok1_temp(end)]
    y_ok1_coef = polyfit(x_pts,y_pts,3);
    y_ok1_int = polyval(y_ok1_coef,xq);
    
    %- Isolate segment
    [c imin] = min(y_ok1_int);
    y_temp = y_ok1_int(1:imin);
    
    dy = diff(y_ok1_temp);
    [c imin] = min(dy);
    boundary_max = y_ok1_temp(imin);
    boundary_min = y_ok1_temp(imin+1);
    
    segment_model = y_temp(y_temp>boundary_min & y_temp<boundary_max);
    x = 1:length(segment_model);
    x1 = x_ok1_temp(imin);x2 = x_ok1_temp(imin+1);
    l = x2-x1-1;
    xq = x(1):(x(end)-x(1))/l :x(end);
    segment_interpolate = interp1(x,segment_model,xq)
    
    y_up = [y_ok1_temp(1:imin) ; segment_interpolate' ; y_ok1_temp(imin+1:end)]
    x_up = [x_ok1_temp(1:imin-1) ; [x1:x2]' ;x_ok1_temp(imin+1:end)]
        
    %- Check First point of upper boundary is on reference line
    if x_up(1)~=x1
        x_diff = x1 - x_up(1);
        y_diff = max(y_up) - y_up(1);
        x_temp = [x1 x1-round(x_diff/2)];
        y_temp = [max(y_up) max(y_up)-y_diff/2];
        x_up = [x_temp';x_up];
        y_up = [y_temp';y_up];
    end
    
function [centreline] = find_centreline(x_ok1,y_ok1,x_ok2,y_ok2,x_centre,z_ref_matrix,flag_plot)
    i=1;
    y_ok1 = -y_ok1+z_ref_matrix; y_ok2 = -y_ok2+z_ref_matrix;
    ind1 = find(y_ok1<0)
    y_ok1(ind1) = [];x_ok1(ind1) = [];
    ind2 = find(y_ok2<0)
    y_ok2(ind2) = [];x_ok2(ind2) = [];
    
    x_ok1 = [x_ok1(1);x_ok1;x_ok1(end)];
    y_ok1 = [0;y_ok1;0];
    x_ok2 = [x_ok2(1);x_ok2;x_ok2(end)];
    y_ok2 = [0;y_ok2;0];
    
    
%     figure(),hold all, plot(x_ok1,y_ok1), plot(x_ok2,y_ok2)
    
    aorta_out = [x_ok1'; y_ok1'];
    aorta_in = [x_ok2'; y_ok2'];
    if flag_plot==1
        figure()
    end
    for k=0:2:180
        if k==90
            ymax = 100;
            line = -9:60;
            l = length(line);
            line_prosp = [x_centre*ones(1,l);line];
        else
            ymax = tan(k*pi/180)*x_centre;
            line = ymax - (ymax/x_centre)*[1:2*x_centre];
            line_prosp = [[1:2*x_centre];line];
        end
        if flag_plot==1
            hold all, plot(x_ok1,y_ok1,'b'),plot(x_ok2,y_ok2,'b'),plot(line_prosp(1,:),line_prosp(2,:));
            ylim([-2 60]);
        end
        Pin = InterX(aorta_out,line_prosp);
        Pout = InterX(aorta_in,line_prosp);
        lin = size(Pin,2); lout = size(Pout,2);
        if lin ~= lout
            if lin>lout
                Pin=[Pin(:,1) Pin(:,end)];
            else
                Pout(:,lout)=[];
            end
        end
        
        if k==0
            %- Point #1
            x_mid = Pout(1,1)+(Pin(1,1)-Pout(1,1))/2;
            y_mid = Pout(2,1)+(Pin(2,1)-Pout(2,1))/2;
            centreline(i,:) = [(x_mid) (y_mid)];
            if flag_plot==1
                plot(centreline(i,1),centreline(i,2),'*k');
            end
            %- Point #2
            x_mid = Pout(1,2)-(Pout(1,2)-Pin(1,2))/2;
            y_mid = Pout(2,2)-(Pout(2,2)-Pin(2,2))/2;
            centreline(i+1,:) = [(x_mid) (y_mid)];
            if flag_plot==1
                plot(centreline(i+1,1),centreline(i+1,2),'*k');
            end
            i = i+2;
        else
            x_mid = Pout(1,1)+(Pin(1,1)-Pout(1,1))/2;
            y_mid = Pout(2,1)+(Pin(2,1)-Pout(2,1))/2;
            centreline(i,:) = [(x_mid) (y_mid)];
            if flag_plot==1
                plot(centreline(i,1),centreline(i,2),'*k');
            end
            i=i+1;
        end
    
%         if size(Pout,2)==1
%             x_mid = Pout(1,1)+(Pin(1,1)-Pout(1,1))/2;
%             y_mid = Pout(2,1)+(Pin(2,1)-Pout(2,1))/2;
%             centreline(i,:) = [(x_mid) (y_mid)];
%             if flag_plot==1
%                 plot(centreline(i,1),centreline(i,2),'*k');
%             end
%             i=i+1;
%         else
%             %- Point #1
%             x_mid = Pout(1,1)+(Pin(1,1)-Pout(1,1))/2;
%             y_mid = Pout(2,1)+(Pin(2,1)-Pout(2,1))/2;
%             centreline(i,:) = [(x_mid) (y_mid)];
%             if flag_plot==1
%                 plot(centreline(i,1),centreline(i,2),'*k');
%             end
%             %- Point #2
%             x_mid = Pout(1,2)-(Pout(1,2)-Pin(1,2))/2;
%             y_mid = Pout(2,2)-(Pout(2,2)-Pin(2,2))/2;
%             centreline(i+1,:) = [(x_mid) (y_mid)];
%             if flag_plot==1
%                 plot(centreline(i+1,1),centreline(i+1,2),'*k');
%             end
%             i = i+2;
%         end    function mouseMove (object, eventdata)
C = get (gca, 'CurrentPoint');
title(gca, ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);
        if flag_plot==1
            plot(Pout(1,:),Pout(2,:),'*r');
            plot(Pin(1,:),Pin(2,:),'*r');
        end
        clear Pout Pin
    end
    
    %% Sort Centreline values
    temp = centreline(2,:);
    centreline(2,:) = [];
    centreline(end+1,:) = temp;
    centreline = round(centreline);

function P = InterX(L1,varargin)
    %...Argument checks and assignment of L2
    error(nargchk(1,2,nargin));
    if nargin == 1,
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

    %...Obtain the segments where an intersectihttps://www.facebook.com/on is expected
    [i,j] = find(C1 & C2); 
    if isempty(i),P = zeros(2,0);return; end;
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
              
function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);

function mask = create_mask(Xf,x,y,line_nb)
    mask = zeros(size(Xf,1),size(Xf,2));
    for k=1:length(x)
        mask(y(k),x(k)) = 1;
    end
    mask(line_nb,:)=1;
    
function dataout = scaledata(datain,minval,maxval)
    dataout = datain - min(datain(:));
    dataout = (dataout/range(dataout(:)))*(maxval-minval);
    dataout = dataout + minval;

function l = calculate_length(x,y,info)
    l_datapoints = length(x);
    for k=2:l_datapoints
        x_temp = x(k)-x(k-1);
        y_temp = y(k)-y(k-1);
        l_temp(k) = sqrt(y_temp^2 + x_temp^2);
    end
    l = sum(l_temp);
    l = l*info.PixelSpacing(1);

function starttrack(imagefig,varargins)
% global axsize
% imagefig.Units = 'pixels';
% imagefig.Position = [imagefig.Position(1) axsize(2) axsize(3) axsize(4)]
% get(0,'PointerLocation')
disp('tracking started')
set(gcf,'windowbuttondownfcn',{@stoptrack},'userdata',[]);
set(gcf,'windowbuttonmotionfcn',...
             'set(gcf,''userdata'',[get(gcf,''userdata'');get(0,''pointerlocation'')])');

function stoptrack(imagefig,varargins)

set(gcf,'windowbuttonmotionfcn',[]);

disp('tracking stopped')

units0=get(0,'units');
unitsf=get(gcf,'units');
unitsa=get(gca,'units');
                           
set(0,'units','pixels');
set(gcf,'units','pixels');
set(gca,'units','pixels');

x=get(gca,'xlim');
y=get(gca,'ylim');
axsize = plotboxpos(gca);
figsize=get(gcf,'position'); 
ratio=diag([diff(x)/axsize(3) diff(y)/axsize(4)]);
shift=figsize(1:2)+[axsize(1) axsize(2)];

set(0,'units',units0);
set(gcf,'units',unitsf);
set(gca,'units',unitsa);

data = get(gcf,'userdata');
mousetrail=(data-repmat(shift,size(data,1),1))*ratio;

set(gcf,'windowbuttondownfcn',{@starttrack});
set(0,'userdata',mousetrail);

function TT = TTAlgorithm(signal,f,algorithm,waveform,continuous,show)
% Description:
%   TTAlgorithm applies user defined pulse wave analysis algorithms in
%   order to determine the transit time between two wave forms.  This
%   software assumes that the two waveforms are measured simultaneously,
%   and no temporal alignment is needed.
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
% ************
% ** Please cite the following paper if used for publication **
% 
% N R Gaddum, J Alastruey, P Beerbaum, P Chowienczyk, T Schaeffter,
% A technical assessment of pulse wave velocity algorithms applied 
% to non-invasive arterial waveforms, Annals of Biomedical Engineering,
% 2013, 41(12):2617-29. DOI: 10.1007/s10439-013-0854-y
% 
% ************
% 
% Author:
% Dr. Nicholas Gaddum
% 
% Revisions:
% 2012-Aug-06   Created function
% 2012-Aug-09   Removed replication of time parameters.  Only frequency 
% 	remaining.  Added pressure/area vs. velocity/flow parameter for
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

%%%%%%% NOTE, wave feet can be poorly located if your frequency is
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
[indmaxs,indmins,gradins] = Find_Landmarks(t,signal,t_intv,npoly_1,npoly1_1,f,algorithm,waveform,continuous,show);

switch algorithm
    case 3
        % Least Squares difference of the systolic upstroke
        t_int = PWV_Calulator_Least_Squares(t,signal,indmaxs,indmins);
    case 4
        % Cross correlation of the entire waveform
        t_int = PWV_Calulator_CC_cycle(t,signal,indmaxs,indmins,gradins);
    otherwise
        % Both foot to foot algorithms
        t_int = PWV_Calulator_FTF(t,signal,indmaxs,indmins,gradins,show);
end
TT = t_int;

function [t_int,t_foot,horint] = PWV_Calulator_FTF(t,signal,indmaxs,indmins,gradins,show)

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
        plot(t_foot,horint,'rx','MarkerSize',12,'Linewidth',2)
end

non_zeros = find(t_foot(1,:).*t_foot(2,:));

t_int = TT;

function [indmaxs,indmins,gradins] = Find_Landmarks(t,signal,t_intv,npoly_1,npoly1_1,f,algorithm,waveform,continuous,show)
%%%%%%%%%%%%% Locate Minima, Maxima and Max Gradients %%%%%%%%%%%%%
t = t';
signal = signal';

switch show
    case 0
    case 1
        fig1 = figure(1);
        set(fig1,'Position',[50 50  826 274])
        set(gca,'Position',[0.06 0.16 0.93 0.72])
        cla
        hold on
        plot(t,signal(:,1),'k')
        plot(t,signal(:,2),'Color',[0.5 0.5 0.5])
        low = min(min(signal));
        high = max(max(signal));
        axis([ t(1) t(end) (low - 0.1*abs(low)) (high + 0.1*abs(high)) ])
        title('TT Analysis of Waveform Data','FontSize',16)
        xlabel('Time (seconds)','FontSize',14)
        ylabel('Waveform Magnitude','FontSize',14)
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
        threshold = 2.5.*sqrt(var( signal(t_intv(2) : t_intv(1) ) ));
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
        npoly1 = npoly1_1;
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
            figure(1)
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

name
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
        figure(1)
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

function pos = plotboxpos(h)
%PLOTBOXPOS Returns the position of the plotted axis region
%
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%
%   h:      axis handle of a 2D axis (if ommitted, current axis is used).
%
% Output variables:
%
%   pos:    four-element position vector, in same units as h

% Copyright 2010 Kelly Kearney

% Check input

if nargin < 1
    h = gca;
end

if ~ishandle(h) || ~strcmp(get(h,'type'), 'axes')
    error('Input must be an axis handle');
end

% Get position of axis in pixels

currunit = get(h, 'units');
axisPos  = getpixelposition(h);

% Calculate box position based axis limits and aspect ratios

darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');

if ~darismanual && ~pbarismanual
    
    pos = axisPos;
    
else

    xlim = get(h, 'XLim');
    ylim = get(h, 'YLim');
    
    % Deal with axis asclimits auto-set via Inf/-Inf use
    
    if any(isinf([xlim ylim]))
        hc = get(h, 'Children');
        hc(~arrayfun( @(h) isprop(h, 'XData' ) & isprop(h, 'YData' ), hc)) = [];
        xdata = get(hc, 'XData');
        if iscell(xdata)
            xdata = cellfun(@(x) x(:), xdata, 'uni', 0);
            xdata = cat(1, xdata{:});
        end
        ydata = get(hc, 'YData');
        if iscell(ydata)
            ydata = cellfun(@(x) x(:), ydata, 'uni', 0);
            ydata = cat(1, ydata{:});
        end
        isplotted = ~isinf(xdata) & ~isnan(xdata) & ...
                    ~isinf(ydata) & ~isnan(ydata);
        xdata = xdata(isplotted);
        ydata = ydata(isplotted);
        if isempty(xdata)
            xdata = [0 1];
        end
        if isempty(ydata)
            ydata = [0 1];
        end
        if isinf(xlim(1))
            xlim(1) = min(xdata);
        end
        if isinf(xlim(2))
            xlim(2) = max(xdata);
        end
        if isinf(ylim(1))
            ylim(1) = min(ydata);
        end
        if isinf(ylim(2))
            ylim(2) = max(ydata);
        end
    end

    dx = diff(xlim);
    dy = diff(ylim);
    dar = get(h, 'DataAspectRatio');
    pbar = get(h, 'PlotBoxAspectRatio');

    limDarRatio = (dx/dar(1))/(dy/dar(2));
    pbarRatio = pbar(1)/pbar(2);
    axisRatio = axisPos(3)/axisPos(4);

    if darismanual
        if limDarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/limDarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * limDarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    elseif pbarismanual
        if pbarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/pbarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * pbarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    end
end

% Convert plot box position to the units used by the axis

hparent = get(h, 'parent');
hfig = ancestor(hparent, 'figure'); % in case in panel or similar
currax = get(hfig, 'currentaxes');

temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', hparent);
set(temp, 'Units', currunit);
pos = get(temp, 'position');
delete(temp);

set(hfig, 'currentaxes', currax);

function [x_left y_left x_right y_right] = coordinate_boundaries_lower(xline, yline, z_ref_matrix, z_ref_diaphragm_matrix)
    [c ind_x] = min(xline)
    xline = [xline(ind_x:end) ; xline(1:ind_x-1)]
    yline = [yline(ind_x:end) ; yline(1:ind_x-1)]
    
    xmin = min(xline); xmax = max(xline); spread = xmax - xmin +1;
    i = 1;
    for x=xmin:xmax
        x_ind = find(xline==x)
        if length(x_ind)==1
            a=1;
        elseif length(x_ind)==2
            x_left(i) = x; y_left(i) = yline(min(x_ind));
            x_right(i) = x; y_right(i) = yline(max(x_ind));
            i = i+1;
        else
            x_left(i) = x; x_right(i) = x;
            y_possible = yline(x_ind);
            [idx C] = kmeans(y_possible,2);
            y_left(i) = round(max(C)); y_right(i) = round(min(C));
            i = i+1;
        end
    end
    
    dy_right = abs(diff(y_right));  ind_right = find(dy_right>5)
    dy_left = abs(diff(y_left));  ind_left = find(dy_left>5)
    ind = [ind_right ; ind_left]
    if ~isempty(ind)
        ind = ind+1;
        x_left(ind) = []; y_left(ind) = [];
        x_right(ind) = []; y_right(ind) = [];
    end
    x1 = x_left; x2 = x_right; y1 = y_left; y2 = y_right;
    y_left = x2; x_left = y2;
    y_right = x1; x_right = y1;
    %- interpolate top of the vessel
    if y_left(1)~=z_ref_matrix
        y_left_temp = [z_ref_matrix:y_left(1)-1]; y_right_temp = [z_ref_matrix:y_right(1)-1];
        x_left_temp = x_left(1)*ones(1,length(y_left_temp)); x_right_temp = x_right(1)*ones(1,length(y_right_temp)); 
        y_left = [y_left_temp y_left]; x_left = [x_left_temp x_left];
        y_right = [y_right_temp y_right]; x_right = [x_right_temp x_right];
    end
    %- interpolate bottom of the vessel
    if y_left(end)~=z_ref_diaphragm_matrix
        y_left_temp = [y_left(end):z_ref_diaphragm_matrix]; y_right_temp = [y_right(end):z_ref_diaphragm_matrix];
        x_left_temp = x_left(end)*ones(1,length(y_left_temp)); x_right_temp = x_right(end)*ones(1,length(y_right_temp));
        y_left = [y_left y_left_temp]; x_left = [x_left x_left_temp];
        y_right = [y_right y_right_temp]; x_right = [x_right x_right_temp];
    end
    a=1;
    
function [Tasc Tdesc Tdia Qasc Qdesc Qdia] = extract_Q_wave(filepath)

%- Load file 
% data = readtable(filepath,'Delimiter','\t','Format','%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
% datacell = data{:,:};
fid = fopen(filepath,'r','n')
bytes = fread(fid)';
asciibytes = bytes(1:2:end) % strip out the zero bytes 
fid = fopen('Temporary.txt','w','n','UTF-8');
fwrite(fid,asciibytes);
fclose(fid);
data = readtable('Temporary.txt','Delimiter','\t','Format','%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
datacell = data{:,:};

%- Retrieve first column
for i=1:size(datacell,1)
    for j=1:size(datacell,2)
        if ~isnan(str2double(datacell{i,j}))
            datacell{i,j} = str2double(datacell{i,j});
        end
        if j==1
            C1{i,1} = datacell{i,1};
        end
    end
end

%- Retrieve row of number
i=1;
for k=1:length(C1)
    if isnumeric(C1{k})==1
        id_ok(i) = k;
        i = i+1;
    end
end

%- Isolate series
d_id_ok = diff(id_ok)
[c ind_separate] = findpeaks(d_id_ok)
if size(ind_separate,1)==1
    ind_separate = ind_separate';
end
ind_separate = [0 ; ind_separate]
for k=1:3
    Time_ind{k} = [id_ok(ind_separate(k)+1):id_ok(ind_separate(k+1))]';
end


%% Retrieve flow waves
idx = find(strcmp(C1,'Time(s)'));
if length(idx)<4
    variable_output = [{'Qasc'},{'Qdesc'}];
    variable_time_output = [{'Tasc'},{'Tdesc'}];
else
    variable_output = [{'Qasc'},{'Qdesc'},{'Qdia'}];
    variable_time_output = [{'Tasc'},{'Tdesc'},{'Tdia'}];
end


for k=1:length(idx)-1
    variables_name = [datacell(idx(1),:)];
    %- Time vector
    eval([variable_time_output{k} '=datacell(Time_ind{k},1);'])
    eval([variable_time_output{k} '= cell2mat(' variable_time_output{k} ')'])
    %- Q vector
    id_var = find(strcmp(variables_name,'Flow Rate(ml/s)'));
    eval([variable_output{k} '=datacell(Time_ind{k},id_var);'])
    eval(['id_ko = cellfun(@ischar,' variable_output{k} ')'])
    eval([variable_output{k} '(id_ko)={0}'])
    eval([variable_output{k} '= cell2mat(' variable_output{k} ')'])
end

if ~exist('Qdia')
    Tdia = []; Qdia = [];
end

function [centerline] = find_centreline_diaphragm(x_left,y_left,x_right,y_right,flag_plot)
    for k=1:length(y_left)
        x(k) = round(mean([x_left(k) x_right(k)]))
        y(k) = y_left(k);
    end
    centerline = [x ; y];
    
    if flag_plot ==1
        figure(), hold all, plot(x_left,y_left,'b'), plot(x_right,y_right,'b'), plot(x,y,'r')
    end
