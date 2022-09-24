function varargout = GUI_coupling_LA_LV(varargin)
% GUI_COUPLING_LA_LV MATLAB code for GUI_coupling_LA_LV.fig
%      GUI_COUPLING_LA_LV, by itself, creates a new GUI_COUPLING_LA_LV or raises the existing
%      singleton*.
%
%      H = GUI_COUPLING_LA_LV returns the handle to a new GUI_COUPLING_LA_LV or the handle to
%      the existing singleton*.
%
%      GUI_COUPLING_LA_LV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_COUPLING_LA_LV.M with the given input arguments.
%
%      GUI_COUPLING_LA_LV('Property','Value',...) creates a new GUI_COUPLING_LA_LV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_coupling_LA_LV_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed t, cla reset;o GUI_coupling_LA_LV_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)". 
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_coupling_LA_LV

% Last Modified by GUIDE v2.5 03-May-2022 15:35:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_coupling_LA_LV_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_coupling_LA_LV_OutputFcn, ...
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


% --- Executes just before GUI_coupling_LA_LV is made visible.
function GUI_coupling_LA_LV_OpeningFcn(hObject, eventdata, handles, varargin)
global a1_color a2_color a3_color a4_color lw handles

% Choose default command line output for GUI_coupling_LA_LV    
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%- Load dataset
% load('data.mat')

%- Plots parameters
a1_color = [100 229 255]/256;
a2_color = [255 243 100]/256;
a3_color = [100 255 101]/256;
a4_color = [255 100 100]/256;
lw = 2;

% UIWAIT makes GUI_coupling_LA_LV wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_coupling_LA_LV_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_compute.
function button_compute_Callback(hObject, eventdata, handles)
global la lv i_LA v_LA  handles Global_Function Reservoir Early_Diastole Conduit Booster i_LA

% %- Clear plots and metrics
% axes(handles.axes_dvLA), cla reset;
% axes(handles.axes_vLA), cla reset;
% axes(handles.axes_vLV), cla reset;
% axes(handles.axes_dvLV), cla reset;

%- Realign signal
[c i] = min(la)
if i>round(length(la)/2)
    la = la(1:i);
    lv = lv(1:i);
end

%- Find landmarks & Metrics
[i_LA v_LA] = LA_landmarks(la);
for i=1:length(i_LA)-1
    if i_LA(i+1) == i_LA(i)
        i_LA(i+1) = i_LA(i+1)+1;
    end
end
refresh()

% --- Executes on button press in button_systole_end.
function button_systole_end_Callback(hObject, eventdata, handles)
global a1_color a2_color a3_color a4_color lw i_LA v_LA la lv Global_Function Reservoir Early_Diastole Conduit Booster i_LA

%- Select point
[xi,yi] = getpts(handles.axes_vLA);
i_LA(1) = round(xi); v_LA(1) = la(i_LA(1));

%- Compute metrics and refresh
refresh()

% --- Executes on button press in button_proto_diastole.
function button_proto_diastole_Callback(hObject, eventdata, handles)
global i_LA v_LA la lv Global_Function Reservoir Early_Diastole Conduit Booster i_LA

%- Select point
[xi,yi] = getpts(handles.axes_vLA);
i_LA(2) = round(xi); v_LA(2) = la(i_LA(2))

%- Compute metrics and refresh
refresh()

% --- Executes on button press in button_mid_systole.
function button_mid_systole_Callback(hObject, eventdata, handles)
global i_LA v_LA la lv Global_Function Reservoir Early_Diastole Conduit Booster i_LA

%- Select point
[xi,yi] = getpts(handles.axes_vLA);
i_LA(3) = round(xi); v_LA(3) = la(i_LA(3))

%- Refresh
refresh()

% --- Executes on button press in button_export.
function button_export_Callback(hObject, eventdata, handles)
global la lv path filename Global_Function Reservoir Early_Diastole Conduit Booster HR i_LA Patient
[file_output,path_output,indx] = uiputfile('.xls');
generate_spreadsheet_LA_LV(file_output,path_output,Patient,Global_Function,Reservoir,Early_Diastole,Conduit,Booster,HR,la,lv,i_LA)
a = 1;

% --- Executes on button press in button_import.
function button_import_Callback(hObject, eventdata, handles)
global la lv path filename HR Patient
% hObject    handle to button_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dname_folder = pwd;
[filename,path]=uigetfile({'*.*','All Files'},'Select flow file for calculating transit time');
[la,lv,HR,Patient] = extract_data_LA_LV(path, filename)
cd(dname_folder)

%% ADDITIONAL FUNCTIONS %%
function [] = refresh()
global handles a1_color a2_color a3_color a4_color lw i_LA v_LA la lv HR Global_Function Reservoir Early_Diastole Conduit Booster

%- Clear plots and metrics
axes(handles.axes_dvLA),
cla reset;
axes(handles.axes_vLA),
cla reset;
axes(handles.axes_vLV),
cla reset;
axes(handles.axes_dvLV),
cla reset;

%- Plot dvLA
axes(handles.axes_dvLA)
hold all,
plot(diff(la),'r','linewidth',lw)
line([0 length(la)],[0 0],'color','k')
xlim([0 length(la)]), yl = ylim,
line([i_LA(1) i_LA(1)],[yl(1) yl(2)],'color','k')
line([i_LA(2) i_LA(2)],[yl(1) yl(2)],'color','k')
line([i_LA(3) i_LA(3)],[yl(1) yl(2)],'color','k')
xt = round(i_LA/length(la)*60/HR*1000)/1000;
set(gca,'XTick',[i_LA])
set(gca,'XTickLabel',xt)
ylabel('dV_L_A volume  (ml.s^-^1)')

%- Plot vLA
axes(handles.axes_vLA)
hold all,
a1 = area(1:i_LA(1),la(1:i_LA(1))), a1.FaceColor = a1_color
a2 = area(i_LA(1):i_LA(2),la(i_LA(1):i_LA(2))), a2.FaceColor = a2_color
a3 = area(i_LA(2):i_LA(3),la(i_LA(2):i_LA(3))), a3.FaceColor = a3_color
a4 = area(i_LA(3):i_LA(4),la(i_LA(3):i_LA(4))), a4.FaceColor = a4_color
plot(la,'r','linewidth',lw),
ylabel('LA volume (ml)'), xlim([0 length(la)]), ylim([min(la)-1 (round(max(la)/10)+1)*10])
set(gca,'XTick',[i_LA])
set(gca,'XTickLabel',xt)

%- Plot vLV
axes(handles.axes_vLV)
hold all,
a1 = area(1:i_LA(1),lv(1:i_LA(1))), a1.FaceColor = a1_color
a2 = area(i_LA(1):i_LA(2),lv(i_LA(1):i_LA(2))), a2.FaceColor = a2_color
a3 = area(i_LA(2):i_LA(3),lv(i_LA(2):i_LA(3))), a3.FaceColor = a3_color
a4 = area(i_LA(3):i_LA(4),lv(i_LA(3):i_LA(4))), a4.FaceColor = a4_color
plot(lv,'r','linewidth',lw),
ylabel('LV volume (ml)'), xlim([0 length(lv)]), ylim([min(lv)-1 (round(max(lv)/10)+1)*10])
set(gca,'XTick',[i_LA])
set(gca,'XTickLabel',xt)

%- Plot dvLV
axes(handles.axes_dvLV)
hold all,
plot(diff(lv),'r','linewidth',lw),line([0 length(lv)],[0 0],'color','k')
xlim([0 length(lv)]), yl = ylim,
line([i_LA(1) i_LA(1)],[yl(1) yl(2)],'color','k')
line([i_LA(2) i_LA(2)],[yl(1) yl(2)],'color','k')
line([i_LA(3) i_LA(3)],[yl(1) yl(2)],'color','k')
set(gca,'XTick',[i_LA])
set(gca,'XTickLabel',xt)
ylabel('dV_L_V volume (ml.s^-^1)')
xlabel('t (s)')

%- Metrics
[Global_Function Reservoir Early_Diastole Conduit Booster] = metrics_LA(v_LA);
set(handles.field_Global_Function, 'String', num2str(Global_Function,2));
set(handles.field_Reservoir, 'String', num2str(Reservoir,3));
set(handles.field_Early_Diastole, 'String', num2str(Early_Diastole,2));
set(handles.field_Conduit, 'String', num2str(Conduit,2));
set(handles.field_Booster, 'String', num2str(Booster,2));

%- Edit Timing
% string1 = sprintf('M(1) = %.3f', M(1));
set(handles.edit_systole, 'String', xt(1));
set(handles.edit_proto_diastole, 'String', xt(2));
set(handles.edit_mid_diastole, 'String', xt(3));

function edit_systole_Callback(hObject, eventdata, handles)
global handles a1_color a2_color a3_color a4_color lw i_LA v_LA la lv HR Global_Function Reservoir Early_Diastole Conduit Booster
val = str2double(get(handles.edit_systole, 'String'));
i_LA(1) = round(val/(60/HR/size(la,1)));
v_LA(1) = la(i_LA(1));
refresh()

% --- Executes during object creation, after setting all properties.
function edit_systole_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_systole (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_proto_diastole_Callback(hObject, eventdata, handles)
global handles a1_color a2_color a3_color a4_color lw i_LA v_LA la lv HR Global_Function Reservoir Early_Diastole Conduit Booster
val = str2double(get(handles.edit_proto_diastole, 'String'));
i_LA(2) = round(val/(60/HR/size(la,1)));
v_LA(2) = la(i_LA(2));
refresh()

% --- Executes during object creation, after setting all properties.
function edit_proto_diastole_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_proto_diastole (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_mid_diastole_Callback(hObject, eventdata, handles)
global handles a1_color a2_color a3_color a4_color lw i_LA v_LA la lv HR Global_Function Reservoir Early_Diastole Conduit Booster
val = str2double(get(handles.edit_mid_diastole, 'String'))
i_LA(3) = round(val/(60/HR/size(la,1)));
v_LA(3) = la(i_LA(3));
refresh()

% --- Executes during object creation, after setting all properties.
function edit_mid_diastole_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mid_diastole (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [LA_vol,LV_vol,HR,Patient] = extract_data_LA_LV(path, filename)

filepath = [path, filename]
fid = fopen(filepath,'r','n');
bytes = fread(fid)';
asciibytes = bytes(1:2:end); % strip out the zero bytes 
fid = fopen('Temporary.txt','w','n','UTF-8');
fwrite(fid,asciibytes);
fclose(fid);
data = readtable('Temporary.txt','Delimiter','\t','Format','%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
datacell = data{:,:};

%% Retrieve first column
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

%% Find 2CV/4CV biplan part
for i=1:length(C1)
    if ~isempty(strfind(C1{i},'2CV/4CV biplan'))
        ind_start_str = i;
    end
end
  
%% Retrieve row of number with useful data
i=1;
for j=1:length(C1)
    if isnumeric(C1{j})==1
       id_ok(i) = j;
       i = i+1;
    end
end
d_id_ok = diff(id_ok); 
[c ind_peak] = findpeaks(d_id_ok); 
index{1} = id_ok(1:ind_peak(1))
for i=1:length(ind_peak)-1
    index{i+1} = id_ok(ind_peak(i)+1:ind_peak(i+1));
end
index{length(ind_peak)+1} = id_ok(ind_peak(end)+1:end);
  
for i=1:length(index)
    id_start_data(i) = index{i}(1);
end

%% HR
i=1;
for j=1:length(C1)
    if ~isempty(strmatch(C1{j},'HR','exact'))
       id_HR(i) = j;
       i = i+1;
    end
end

%% Patient info
i=1;
for j=1:length(C1)
    if ~isempty(strmatch(C1{j},'Patient','exact'))
       id_Patient = j;
       Patient.Name = datacell{id_Patient,2}
       i = i+1;
    end
end

i=1;
for j=1:length(C1)
    if ~isempty(strmatch(C1{j},'Birth Date','exact'))
       id_DOB = j;
       Patient.DOB = datacell{id_DOB,2}
       i = i+1;
    end
end


i=1;
for j=1:length(C1)
    if ~isempty(strmatch(C1{j},'PatientID','exact'))
       id_PatientID = j;
       Patient.PatientID = datacell{id_PatientID,2}
       i = i+1;
    end
end

i=1;
for j=1:length(C1)
    if ~isempty(strmatch(C1{j},'Study Date','exact'))
       id_ScanDate = j;
       Patient.ScanDate = datacell{id_ScanDate,2}
       i = i+1;
    end
end

i=1;
for j=1:length(C1)
    if ~isempty(strmatch(C1{j},'AccessionNumber','exact'))
       id_AccessionNumber(i) = j;
       Patient.AccessionNumber = datacell{id_AccessionNumber,2}
       i = i+1;
    end
end


%% Identify the right row part corresponding to 4CV biplan
[c ind_right] = min(abs(id_start_data-ind_start_str));
id_ok = index{ind_right};

%- Find right HR
[c i] = min(abs(id_HR - id_ok(1)))
id_HR = id_HR(i)
 
%- LA and LV string labels
l=1;
for j=1:size(datacell,2)
    if ~isempty(strfind(datacell{id_ok(1)-2,j},'la'))
       ind_label_la(l) = j;
       l=l+1;
    end
end

l=1
for j=1:size(datacell,2)
    if ~isempty(strmatch(datacell{id_ok(1)-2,j},'lv','exact'))
       ind_label_lv(l) = j;
       l=l+1;
    end
end

%- Decide which LA and LV columns to use from the spreadsheet
ind_label_la = ind_label_la(ind_label_la>ind_label_lv(end))
ind_label = [ind_label_lv(end) ind_label_la(1) ]

%- Retrieve all LA/LV data
for i=1:length(ind_label)
    m=1;
    for j=1:length(id_ok)
        LA_data(m,i) = datacell{id_ok(j),ind_label(i)};
        m=m+1;
    end
end
LV_vol = LA_data(:,1);
LA_vol = LA_data(:,2);
HR = datacell{id_HR,2}

function [i v] = LA_landmarks(la)
    
    %- Peak volume
    [v(1) i(1)] = max(la);
    
    %- Start conduit volume
    dla = diff(la);
    [ind_intersection] = intersection(dla);
    i_LA_start_conduit = ind_intersection(ind_intersection>i(1));
    i(2) = i_LA_start_conduit(1);
    v(2) = la(i(2));
    
    %- End conduit volume
    i(3) = ind_intersection(end-1);
    v(3) = la(i(3));
    
    %- End
    i(4) = length(la);
    v(4) = la(end);
    
function [ind_intersection] = intersection(flow)
    l = length(flow);
    i = 1;
    for k=2:l
        test = flow(k-1)*flow(k);
        if test<0
            ind_intersection(i) = k;
            i = i+1;
        end
    end
    
function [Global_Function Reservoir Early_Diastole Conduit Booster] = metrics_LA(v_LA)
    Global_Function = (v_LA(1)-v_LA(4))*100/v_LA(1);
    Reservoir = (v_LA(1)-v_LA(4))*100/v_LA(4);
    Early_Diastole = (v_LA(1)-v_LA(2))*100/v_LA(1);
    Conduit = (v_LA(1)-v_LA(3))*100/v_LA(1);
    Booster = (v_LA(3)-v_LA(4))*100/v_LA(3);
    
function [] = generate_spreadsheet_LA_LV(filename,path,Patient,Global_Function,Reservoir,Early_Diastole,Conduit,Booster,HR,la,lv,i_LA)

fill = {'','','','','',''};

%% PATIENT INFO
%- Check fields exist
if isfield(Patient,'Name')
    Field_Name = [Patient.Name];
else
    Field_Name = '';
end
if isfield(Patient,'PatientID')
    Field_PatientID = [Patient.PatientID];
else
    Field_PatientID = '';
end
if isfield(Patient,'ScanDate')
    Field_ScanDate = [Patient.ScanDate];
else
    Field_ScanDate = '';
end
if isfield(Patient,'DOB') 
    Field_DOB = [Patient.DOB];
else
    Field_DOB = '';
end
if isfield(Patient,'AccessionNumber')
    Field_AccessionNumber = [Patient.AccessionNumber];
else
    Field_AccessionNumber = '';
end

line_PATIENT_INFO = [{'%%% PATIENT INFO','','','','',''};
       {'Name:',Field_Name,'','','',''};
       {'Patient ID:',Field_PatientID,'','','',''};
       {'Scan date:',Field_ScanDate,'','','',''};
       {'DOB:',Field_DOB,'','','',''};
       {'Accession number:',Field_AccessionNumber,'','','',''};
       ];

%% LA PERFORMANCE
line_LA_PERFORMANCE = [{'%%% LA PERFORMANCE','','','','',''};
       {'Global function:',round(100*Global_Function)/100,'%','','',''};
       {'Reservoir:',round(100*Reservoir)/100,'%','','',''};
       {'Early Diastole:',round(100*Early_Diastole)/100,'%','','',''};
       {'Conduit:',round(100*Conduit)/100,'%','','',''};
       {'Booster:',round(100*Booster)/100,'%','','',''};
       ];

%% VOLUMEM
line_VOLUME = [{'','Reservoir','Early Diastole','Conduit','Booster',''};
       {'LA:',la(i_LA(1)),la(i_LA(2)),la(i_LA(3)),la(i_LA(4)),'ml'};
       {'LV:',lv(i_LA(1)),lv(i_LA(2)),lv(i_LA(3)),lv(i_LA(4)),'ml'};
       ];
   
%% DERIVATIVE
%- Reservoir
[val_max_dla_res i_max_dla_res] = max(diff(la(1:i_LA(1))));
t_max_dla_res = i_max_dla_res*60/HR/size(lv,1);
[val_min_dla_res i_min_dla_res] = min(diff(la(1:i_LA(1))));
t_min_dla_res = i_min_dla_res*60/HR/size(lv,1);
[val_max_dlv_res i_max_dlv_res] = max(diff(lv(1:i_LA(1))));
t_max_dlv_res = i_max_dlv_res*60/HR/size(lv,1);
[val_min_dlv_res i_min_dlv_res] = min(diff(lv(1:i_LA(1))));
t_min_dlv_res = i_min_dlv_res*60/HR/size(lv,1);

%- Early Diastole
[val_max_dla_ed i_max_dla_ed] = max(diff(la(i_LA(1)+1:i_LA(2))));
t_max_dla_ed = (i_LA(1)+i_max_dla_ed)*60/HR/size(lv,1);
[val_min_dla_ed i_min_dla_ed] = min(diff(la(i_LA(1)+1:i_LA(2))));
t_min_dla_ed = (i_LA(1)+i_min_dla_ed)*60/HR/size(lv,1);
[val_max_dlv_ed i_max_dlv_ed] = max(diff(lv(i_LA(1)+1:i_LA(2))));
t_max_dlv_ed = (i_LA(1)+i_max_dlv_ed)*60/HR/size(lv,1);
[val_min_dlv_ed i_min_dlv_ed] = min(diff(lv(i_LA(1)+1:i_LA(2))));
t_min_dlv_ed = (i_LA(1)+i_min_dlv_ed)*60/HR/size(lv,1);

%- Conduit
[val_max_dla_con i_max_dla_con] = max(diff(la(i_LA(2)+1:i_LA(3))));
t_max_dla_con = (i_LA(2)+i_max_dla_con)*60/HR/size(lv,1);
[val_min_dla_con i_min_dla_con] = min(diff(la(i_LA(2)+1:i_LA(3))));
t_min_dla_con = (i_LA(2)+i_min_dla_con)*60/HR/size(lv,1);
[val_max_dlv_con i_max_dlv_con] = max(diff(lv(i_LA(2)+1:i_LA(3))));
t_max_dlv_con = (i_LA(2)+i_max_dlv_con)*60/HR/size(lv,1);
[val_min_dlv_con i_min_dlv_con] = min(diff(lv(i_LA(2)+1:i_LA(3))));
t_min_dlv_con = (i_LA(2)+i_min_dlv_con)*60/HR/size(lv,1);

%- Booster
[val_max_dla_boo i_max_dla_boo] = max(diff(la(i_LA(3)+1:i_LA(4))));
t_max_dla_boo = (i_LA(3)+i_max_dla_boo)*60/HR/size(lv,1);
[val_min_dla_boo i_min_dla_boo] = min(diff(la(i_LA(3)+1:i_LA(4))));
t_min_dla_boo = (i_LA(3)+i_min_dla_boo)*60/HR/size(lv,1);
[val_max_dlv_boo i_max_dlv_boo] = max(diff(lv(i_LA(3)+1:i_LA(4))));
t_max_dlv_boo = (i_LA(3)+i_max_dlv_boo)*60/HR/size(lv,1);
[val_min_dlv_boo i_min_dlv_boo] = min(diff(lv(i_LA(3)+1:i_LA(4))));
t_min_dlv_boo = (i_LA(3)+i_min_dlv_boo)*60/HR/size(lv,1);

line_DERIVATIVE = [{'Derivative','','','','',''};
       {'Value max dV LA:',val_max_dla_res,val_max_dla_ed,val_max_dla_con,val_max_dla_boo,'ml.s-1'};
       {'Timing max dV LA:',round(100*t_max_dla_res),round(100*t_max_dla_ed),round(100*t_max_dla_con),round(100*t_max_dla_boo),'s'};
       {'Value min dV LA:',val_min_dla_res,val_min_dla_ed,val_min_dla_con,val_min_dla_boo,'ml.s-1'};
       {'Timing min dV LA:',round(100*t_min_dla_res),round(100*t_min_dla_ed),round(100*t_min_dla_con),round(100*t_min_dla_boo),'s'};
       {'Value max dV LV:',val_max_dlv_res,val_max_dlv_ed,val_max_dlv_con,val_max_dlv_boo,'ml.s-1'};
       {'Timing max dV LV:',round(100*t_max_dlv_res),round(100*t_max_dlv_ed),round(100*t_max_dlv_con),round(100*t_max_dlv_boo),'s'};
       {'Value min dV LV:',val_min_dlv_res,val_min_dlv_ed,val_min_dlv_con,val_min_dlv_boo,'ml.s-1'};
       {'Timing min dV LV:',round(100*t_min_dlv_res),round(100*t_min_dlv_ed),round(100*t_min_dlv_con),round(100*t_min_dlv_boo),'s'};
       ];

%% AREA UNDER THE CURVE
t = [1:size(la,1)]/HR;
line_AUC = [{'Area under the curve','','','','',''};
       {'LA:',round(100*trapz(t(1:i_LA(1)),la(1:i_LA(1)))),round(100*trapz(t(i_LA(1):i_LA(2)),la(i_LA(1):i_LA(2)))),round(100*trapz(t(i_LA(2):i_LA(3)),la(i_LA(2):i_LA(3)))),round(100*trapz(t(i_LA(3):i_LA(4)),la(i_LA(3):i_LA(4)))),''};
       {'LV:',round(100*trapz(t(1:i_LA(1)),lv(1:i_LA(1)))),round(100*trapz(t(i_LA(1):i_LA(2)),lv(i_LA(1):i_LA(2)))),round(100*trapz(t(i_LA(2):i_LA(3)),lv(i_LA(2):i_LA(3)))),round(100*trapz(t(i_LA(3):i_LA(4)),lv(i_LA(3):i_LA(4)))),''};
       ];

%% RR INTERVAL TIME
line_RR_INTERVAL = {'RR interval time',round(100*60/HR)/100,'s','','',''};

%% TIME VARIATIONS
time = ([1:size(la,1)]/size(la,1))*60/HR
dla = [0;diff(la)]
dlv = [0;diff(lv)]
filler = strings(length(dla),1)
line_WAVEFORMS = [["Time (s)","Volume LA (ml)","dV LA (ml.s-1)","Volume LV (ml)","dV LV(ml.s-1)",""];[time',la,dla,lv,dlv,filler]]
      

%% ASSEMBLE TABLE
T = [line_PATIENT_INFO ; fill ; line_LA_PERFORMANCE ; fill ; line_VOLUME ; fill ; line_DERIVATIVE ; fill ; line_AUC ; fill ; line_RR_INTERVAL ; fill ; line_WAVEFORMS ];

%% WRITE FILE
writematrix(T,[path filename]);
