function varargout = GUI_pick_coronal_mode(varargin)
% GUI_PICK_CORONAL_MODE MATLAB code for GUI_pick_coronal_mode.fig
%      GUI_PICK_CORONAL_MODE, by itself, creates a new GUI_PICK_CORONAL_MODE or raises the existing
%      singleton*.
%
%      H = GUI_PICK_CORONAL_MODE returns the handle to a new GUI_PICK_CORONAL_MODE or the handle to
%      the existing singleton*.
%
%      GUI_PICK_CORONAL_MODE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PICK_CORONAL_MODE.M with the given input arguments.
%
%      GUI_PICK_CORONAL_MODE('Property','Value',...) creates a new GUI_PICK_CORONAL_MODE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_pick_coronal_mode_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_pick_coronal_mode_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_pick_coronal_mode

% Last Modified by GUIDE v2.5 21-Apr-2021 16:48:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_pick_coronal_mode_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_pick_coronal_mode_OutputFcn, ...
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


% --- Executes just before GUI_pick_coronal_mode is made visible.
function GUI_pick_coronal_mode_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_pick_coronal_mode (see VARARGIN)

% Choose default command line output for GUI_pick_coronal_mode
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_pick_coronal_mode wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_pick_coronal_mode_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in button_load.
function button_load_Callback(hObject, eventdata, handles)
global z_ref z_ref_diaphragm dname_folder Tasc Tdesc Tdia Qasc Qdesc Qdia...
    Aasc Adesc Adia path_flow flag_diapgragm Im_asc_desc_ref ObjectPolarity_mode...
    centers_asc centers_desc radii_asc radii_desc
%% Select Flow files and transit time
    path_length = []; transit_time = []; centers_asc = []; centers_desc = []; radii_asc = []; radii_desc = [];
    dname_folder = pwd;
    [filename1,path_flow]=uigetfile({'*.*','All Files'},'Select flow file for calculating transit time');
    
    %- Find transit timesetappdata(0,'flow',flow)
    [Tasc Tdesc Tdia Qasc Qdesc Qdia Aasc Adesc Adia] = extract_Q_wave([path_flow filename1]);
    cd(path_flow)
    files = dir('*.dcm');
    info_axial = dicominfo(files(1).name);
    setappdata(0,'info',info_axial)
    setappdata(0,'path_flow',path_flow)
    cd(dname_folder)
    if Tasc(1) == 0
        freq = 1e3/Tasc(2);
    else
        freq = 1e3/Tasc(1);
    end
    %- Open GUI for calculating transit time
    if isempty(Tdia)
        flow = [Qasc Qdesc];
        setappdata(0,'flow',flow)
        setappdata(0,'freq',freq)
        GUI_Transit_Time
    else
        flow = [Qasc Qdesc Qdia];
        area = [Aasc Adesc Adia];
        setappdata(0,'flow',flow)
        setappdata(0,'area',area)
        setappdata(0,'freq',freq)
        GUI_Transit_Time_3_signals
    end

    %% Get reference slices informations
    %- Slice representing ascending and descending cross-section
    z_ref = info_axial.SliceLocation;
    cd(path_flow)
    %- Slice representing only diaphragm section
    if isempty(Tdia)~=1
        [filename2,filepath2]=uigetfile({'*.*','All Files'},...
        'Select Diaphragm Flow Folder');
        cd(filepath2)
        files = dir('*.dcm');
        info_axial = dicominfo(files(1).name);
        z_ref_diaphragm = info_axial.SliceLocation;
        flag_diapgragm = 1;
    else
        flag_diapgragm = 0;
    end
    
    %% Display reference ascending-descending slice
    current_folder = cd;
    dname = uigetdir(path_flow,'Select folder containing coronal plane scans');
    cd(dname)
    [V,spatial,dim] = dicomreadVolume(dname);
    V = squeeze(V);
    info = get_info(dname);
    info_global = info{1};
    
    %- Sort out dimensions
    [x y z x0 y0 z0 x_length y_length z_length pix_space Im_axial] = prepare_volume(V,info);
    [c k_ref_asc] = min(abs(z-z_ref));
    [c k_ref_dia] = min(abs(z-z_ref_diaphragm));
    l = k_ref_dia-k_ref_asc+26;
    k_ref_asc_plot = k_ref_dia-k_ref_asc+6;
    col = 7;
    row = find_row(l);
   
    %- Display
    Im_asc_desc_ref = Im_axial{k_ref_asc};
    axes(handles.axes1), imshow(Im_asc_desc_ref)
    ObjectPolarity_mode = decide_ObjectPolarity_mode(Im_asc_desc_ref);
    setappdata(0,'ObjectPolarity_mode',ObjectPolarity_mode)
    setappdata(0,'k_ref_asc_plot',k_ref_asc_plot)
    setappdata(0,'k_ref_asc',k_ref_asc)
    setappdata(0,'k_ref_dia',k_ref_dia)
    setappdata(0,'col',col)
    setappdata(0,'row',row)
    setappdata(0,'info',info)
    setappdata(0,'info_global',info_global)
    setappdata(0,'pix_space',pix_space)
    setappdata(0,'Im_axial',Im_axial)
    setappdata(0,'x',x)
    setappdata(0,'y',y)
    setappdata(0,'z',z)
    setappdata(0,'V',V)
    setappdata(0,'z_ref',z_ref)
    setappdata(0,'z_ref_diaphragm',z_ref_diaphragm)
    setappdata(0,'flag_diapgragm',flag_diapgragm)
    cd(dname_folder)
    %% Update
    guidata(hObject,handles);

% --- Executes on button press in button_automatic_selection.
function button_automatic_selection_Callback(hObject, eventdata, handles)
    global Im_asc_desc_ref centers_desc radii_desc centers_asc radii_asc  centers_desc_ref radii_desc_ref centers_asc_ref radii_asc_ref dname_folder
    cd(dname_folder)
    centers_desc_ref = []; radii_desc_ref = []; centers_asc_ref = []; radii_asc_ref = [];
    [centers_desc radii_desc centers_asc radii_asc] = find_desc_asc_from_ref_scans(Im_asc_desc_ref,centers_desc_ref,radii_desc_ref,centers_asc_ref,radii_asc_ref);
    axes(handles.axes1),imshow(Im_asc_desc_ref), viscircles(centers_desc, radii_desc,'Color','b'), viscircles(centers_asc, radii_asc,'Color','r');
    %- Store variable in handle
    guidata(hObject, handles);

% --- Executes on button press in button_asc.
function button_asc_Callback(hObject, eventdata, handles)
    global Im_asc_desc_ref seg_I centers_desc radii_desc centers_asc radii_asc;
    
    centers_asc = [];
    imshow(Im_asc_desc_ref)
    if ~isempty(centers_desc)
        viscircles(centers_desc, radii_desc,'Color','b');
    end
    set(gcf,'windowbuttondownfcn',{@starttrack});
    waitfor(0,'userdata')
    asc_boundaries=get(0,'userdata');
    x_asc_boundaries = asc_boundaries(:,1); 
    y_asc_boundaries = handles.axes1.YLim(2)-asc_boundaries(:,2); 
    [xo yo radii_asc] = circle_fit(x_asc_boundaries,y_asc_boundaries);
    viscircles([xo yo], radii_asc,'Color','r')
    centers_asc = [xo yo];
 
% --- Executes on button press in button_desc.
function button_desc_Callback(hObject, eventdata, handles)
    global Im_asc_desc_ref seg_I centers_desc radii_desc centers_asc radii_asc;
    
    centers_desc = [];
    imshow(Im_asc_desc_ref)
    if ~isempty(centers_asc)
        viscircles(centers_asc, radii_asc,'Color','r');
    end
    set(gcf,'windowbuttondownfcn',{@starttrack});
    waitfor(0,'userdata')
    desc_boundaries=get(0,'userdata');
    x_desc_boundaries = desc_boundaries(:,1); 
    y_desc_boundaries = handles.axes1.YLim(2)-desc_boundaries(:,2); 
    [xo yo radii_desc] = circle_fit(x_desc_boundaries,y_desc_boundaries);
    viscircles([xo yo], radii_desc,'Color','b')
    centers_desc = [xo yo];
    
% --- Executes on button press in button_automatic_mode.
function button_automatic_mode_Callback(hObject, eventdata, handles)
    global Im_asc_desc_ref seg_I centers_desc radii_desc centers_asc radii_asc;
    outputname = {'Im_asc_desc_ref','centers_desc','radii_desc','centers_asc','radii_asc'};
    for k=1:length(outputname)
        eval(['setappdata(0,''' outputname{k} ''',' outputname{k} ');'])
    end
    GUI_PWV_coronal_auto_choice

% --- Executes on button press in button_manual_mode.
function button_manual_mode_Callback(hObject, eventdata, handles)
    global Im_asc_desc_ref seg_I centers_desc radii_desc centers_asc radii_asc;
    outputname = {'Im_asc_desc_ref','centers_desc','radii_desc','centers_asc','radii_asc'};
    for k=1:length(outputname)
        eval(['setappdata(0,''' outputname{k} ''',' outputname{k} ');'])
    end
    GUI_PWV_coronal_manual_choice

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADDITIONAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tasc Tdesc Tdia Qasc Qdesc Qdia Aasc Adesc Adia] = extract_Q_wave(filepath)

%- Load file 
% data = readtable(filepath,'Delimiter','\t','Format','%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
% datacell = data{:,:};
fid = fopen(filepath,'r','n');
bytes = fread(fid)';
asciibytes = bytes(1:2:end); % strip out the zero bytes 
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
d_id_ok = diff(id_ok);
[c ind_separate] = findpeaks(d_id_ok);
if size(ind_separate,1)==1
    ind_separate = ind_separate';
end
% if length(ind_separate)==3
%     ind_separate = [0 ; ind_separate ; length(d_id_ok)];
%     for k=1:4
%         Time_ind{k} = [id_ok(ind_separate(k)+1):id_ok(ind_separate(k+1))]';
%     end
% else
%     ind_separate = [0 ; ind_separate ; length(d_id_ok)];
%     for k=1:2
%         Time_ind{k} = [id_ok(ind_separate(k)+1):id_ok(ind_separate(k+1))]';
%     end
% end

if length(ind_separate)==1
    ind_separate = [0 ; ind_separate ; length(id_ok)];
    for k=1:2
        Time_ind{k} = [id_ok(ind_separate(k)+1):id_ok(ind_separate(k+1))]';
    end
else
    ind_separate = [0 ; ind_separate ; length(id_ok)];
    for k=1:length(ind_separate)-1
        Time_ind{k} = [id_ok(ind_separate(k)+1):id_ok(ind_separate(k+1))]';
    end
end

%% Retrieve flow waves
idx1 = find(strcmp(C1,{'Time(s)'}));
idx2 = find(strcmp(C1,{'Time(ms)'}));
idx = [idx1 idx2];

if length(idx)==2
        variable_output_A = [{'Aasc'},{'Adesc'}];
        variable_output_Q = [{'Qasc'},{'Qdesc'}];
        variable_output_t = [{'Tasc'},{'Tdesc'}];
        for k=1:2
            variables_name = [datacell(idx(1),:)];
            for i=1:length(variables_name)
                if length(variables_name{i})>0
                    variables_name{i} = variables_name{i}(1:4);
                end
            end
            %- Time vector
            eval([variable_output_t{k} '=datacell(Time_ind{k},1);'])
            eval([variable_output_t{k} '= cell2mat(' variable_output_t{k} ');'])
            %- Q vector
            id_var = find(strcmp(variables_name,'Flow'));
            eval([variable_output_Q{k} '=datacell(Time_ind{k},id_var);'])
            eval(['id_ko = cellfun(@ischar,' variable_output_Q{k} ');'])
            eval([variable_output_Q{k} '(id_ko)={0};'])
            eval([variable_output_Q{k} '= cell2mat(' variable_output_Q{k} ');'])
            %- A vector
            id_var = find(strcmp(variables_name,'Area'));
            eval([variable_output_A{k} '=datacell(Time_ind{k},id_var);'])
            eval(['id_ko = cellfun(@ischar,' variable_output_A{k} ');'])
            eval([variable_output_A{k} '(id_ko)={0};'])
            eval([variable_output_A{k} '= cell2mat(' variable_output_A{k} ');'])
        end
        %- Make sure all vectors have the same length
        if length(Tasc)~=length(Tdesc)
            [c ind] = min([length(Tasc) length(Tdesc)]) 
            if ind == 2
                Aasc(c+1:end)=[];
                Qasc(c+1:end)=[];
                Tasc(c+1:end)=[];
            else
                Adesc(c+1:end)=[];
                Qdesc(c+1:end)=[];
                Tdesc(c+1:end)=[];
            end          
        end
elseif length(idx)==3
        variable_output_A = [{'Aasc'},{'Adesc'},{'Adia'}];
        variable_output_Q = [{'Qasc'},{'Qdesc'},{'Qdia'}];
        variable_output_t = [{'Tasc'},{'Tdesc'},{'Tdia'}];
        for k=1:3
            variables_name = [datacell(idx(1),:)];
            for i=1:length(variables_name)
                if length(variables_name{i})>0
                    variables_name{i} = variables_name{i}(1:4);
                end
            end
            %- Time vector
            eval([variable_output_t{k} '=datacell(Time_ind{k},1);'])
            eval([variable_output_t{k} '= cell2mat(' variable_output_t{k} ');'])
            %- Q vector
            id_var = find(strcmp(variables_name,'Flow'));
            eval([variable_output_Q{k} '=datacell(Time_ind{k},id_var);'])
            eval(['id_ko = cellfun(@ischar,' variable_output_Q{k} ');'])
            eval([variable_output_Q{k} '(id_ko)={0};'])
            eval([variable_output_Q{k} '= cell2mat(' variable_output_Q{k} ');'])
            %- A vector
            id_var = find(strcmp(variables_name,'Area'));
            eval([variable_output_A{k} '=datacell(Time_ind{k},id_var);'])
            eval(['id_ko = cellfun(@ischar,' variable_output_A{k} ');'])
            eval([variable_output_A{k} '(id_ko)={0};'])
            eval([variable_output_A{k} '= cell2mat(' variable_output_A{k} ');'])
        end
        %- Find the right diaphragm vectors
        if Qdia == 0
            Adia = [];
            Qdia = [];
            Tdia = [];
            l = [length(Tasc) length(Tdesc)]
        else
            l = [length(Tasc) length(Tdesc) length(Tdia)]
        end
        %- Make sure all vectors have the same length
        if length(unique(l))~=1
            ind = min(l);
            %- Ascending aorta signals
            Aasc(ind+1:end)=[];
            Qasc(ind+1:end)=[];
            Tasc(ind+1:end)=[];
            %- Descending aorta signals
            Adesc(ind+1:end)=[];
            Qdesc(ind+1:end)=[];
            Tdesc(ind+1:end)=[];
            %- Diaphragm aorta signals
            Adia(ind+1:end)=[];
            Qdia(ind+1:end)=[];
            Tdia(ind+1:end)=[];     
        end
else
        variable_output_A = [{'Aasc'},{'Adesc'},{'Adia1'},{'Adia2'}];
        variable_output_Q = [{'Qasc'},{'Qdesc'},{'Qdia1'},{'Qdia2'}];
        variable_output_t = [{'Tasc'},{'Tdesc'},{'Tdia1'},{'Tdia2'}];
        for k=1:4
            variables_name = [datacell(idx(1),:)];
            for i=1:length(variables_name)
                if length(variables_name{i})>0
                    variables_name{i} = variables_name{i}(1:4);
                end
            end
            %- Time vector
            eval([variable_output_t{k} '=datacell(Time_ind{k},1);'])
            eval([variable_output_t{k} '= cell2mat(' variable_output_t{k} ');'])
            %- Q vector
            id_var = find(strcmp(variables_name,'Flow'));
            eval([variable_output_Q{k} '=datacell(Time_ind{k},id_var);'])
            eval(['id_ko = cellfun(@ischar,' variable_output_Q{k} ');'])
            eval([variable_output_Q{k} '(id_ko)={0};'])
            eval([variable_output_Q{k} '= cell2mat(' variable_output_Q{k} ');'])
            %- A vector
            id_var = find(strcmp(variables_name,'Area'));
            eval([variable_output_A{k} '=datacell(Time_ind{k},id_var);'])
            eval(['id_ko = cellfun(@ischar,' variable_output_A{k} ');'])
            eval([variable_output_A{k} '(id_ko)={0};'])
            eval([variable_output_A{k} '= cell2mat(' variable_output_A{k} ');'])
        end
        %- Find the right diaphragm vectors
        if Qdia1 ~= 0
            Adia = Adia1;
            Qdia = Qdia1;
            Tdia = Tdia1;
        else
            Adia = Adia2;
            Qdia = Qdia2;
            Tdia = Tdia2;
        end
        %- Make sure all vectors have the same length
        l = [length(Tasc) length(Tdesc) length(Tdia)]
        if length(unique(l))~=1
           ind = min(l);
           %- Ascending aorta signals
           Aasc(ind+1:end)=[];
           Qasc(ind+1:end)=[];
           Tasc(ind+1:end)=[];
           %- Descending aorta signals
           Adesc(ind+1:end)=[];
           Qdesc(ind+1:end)=[];
           Tdesc(ind+1:end)=[];
           %- Diaphragm aorta signals
           Adia(ind+1:end)=[];
           Qdia(ind+1:end)=[];
           Tdia(ind+1:end)=[];     
        end
end

if mean(Qasc)<0
    [Qdesc Qasc] = deal(Qasc,Qdesc);
    [Adesc Aasc] = deal(Aasc,Adesc);
    [Tdesc Tasc] = deal(Tasc,Tdesc);
    
    [Qdesc Qdia] = deal(Qdia,Qdesc);
    [Adesc Adia] = deal(Adia,Adesc);
    [Tdesc Tdia] = deal(Tdia,Tdesc);
end

if ~exist('Qdia')
    Tdia = []; Qdia = []; Adia = [];
end

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

%%%%%%%%%%%%% Minima, Maxima and Max Gradients Sorting - Start %%%%%%%%%%%%%%%%%
% Ensure that similar cycles are compared
trash1 = find(gradins(1,:,1));
trash2 = find(gradins(2,:,1));

ind1 = gradins(1,trash1,1);
ind2 = gradins(2,trash2,1);

% Determine the standard deviation between the number of data points between 
% succfunction [indmaxs,indmins,gradins] = Find_Landmarks(t,signal,t_intv,npoly_1,npoly1_1,f,algorithm,waveform,continuous,show)
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
diff1_stfunction [indmaxs,indmins,gradins] = Find_Landmarks(t,signal,t_intv,npoly_1,npoly1_1,f,algorithm,waveform,continuous,show)
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
signal = signal1';d = std(diff1);

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
signal = signal1';essive cycles for both signals.  Which ever has the lowest standard
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

function hfig = tightfig(hfig)
% tightfig: Alters a figure so that it has the minimum size necessary to
% enclose all axes in the figure without excess space around them.
% 
% Note that tightfig will expand the figure to completely encompass all
% axes if necessary. If any 3D axes are present which have been zoomed,
% tightfig will produce an error, as these cannot easily be dealt with.
% 
% hfig - handle to figure, if not supplied, the current figure will be used
% instead.

    if nargin == 0
        hfig = gcf;
    end

    % There can be an issue with tightfig when the user has been modifying
    % the contnts manually, the code below is an attempt to resolve this,
    % but it has not yet been satisfactorily fixed
%     origwindowstyle = get(hfig, 'WindowStyle');
    set(hfig, 'WindowStyle', 'normal');
    
    % 1 point is 0.3528 mm for future use

    % get all the axes handles note this will also fetch legends and
    % colorbars as well
    hax = findall(hfig, 'type', 'axes');
    
    % get the original axes units, so we can change and reset these again
    % later
    origaxunits = get(hax, 'Units');
    
    % change the axes units to cm
    set(hax, 'Units', 'centimeters');
    
    % get various position parameters of the axes
    if numel(hax) > 1
%         fsize = cell2mat(get(hax, 'FontSize'));
        ti = cell2mat(get(hax,'TightInset'));
        pos = cell2mat(get(hax, 'Position'));
    else
%         fsize = get(hax, 'FontSize');
        ti = get(hax,'TightInset');
        pos = get(hax, 'Position');
    end
    
    % ensure very tiny border so outer box always appears
    ti(ti < 0.1) = 0.15;
    
    % we will check if any 3d axes are zoomed, to do this we will check if
    % they are not being viewed in any of the 2d directions
    views2d = [0,90; 0,0; 90,0];
    
    for i = 1:numel(hax)
        
        set(hax(i), 'LooseInset', ti(i,:));
%         set(hax(i), 'LooseInset', [0,0,0,0]);
        
        % get the current viewing angle of the axes
        [az,el] = view(hax(i));
        
        % determine if the axes are zoomed
        iszoomed = strcmp(get(hax(i), 'CameraViewAngleMode'), 'manual');
        
        % test if we are viewing in 2d mode or a 3d view
        is2d = all(bsxfun(@eq, [az,el], views2d), 2);
               
        if iszoomed && ~any(is2d)
           error('TIGHTFIG:haszoomed3d', 'Cannot make figures containing zoomed 3D axes tight.') 
        end
        
    end
    
    % we will move all the axes down and to the left by the amount
    % necessary to just show the bottom and leftmost axes and labels etc.
    moveleft = min(pos(:,1) - ti(:,1));
    
    movedown = min(pos(:,2) - ti(:,2));
    
    % we will also alter the height and width of the figure to just
    % encompass the topmost and rightmost axes and lables
    figwidth = max(pos(:,1) + pos(:,3) + ti(:,3) - moveleft);
    
    figheight = max(pos(:,2) + pos(:,4) + ti(:,4) - movedown);
    
    % move all the axes
    for i = 1:numel(hax)
        
        set(hax(i), 'Position', [pos(i,1:2) - [moveleft,movedown], pos(i,3:4)]);
        
    end
    
    origfigunits = get(hfig, 'Units');
    
    set(hfig, 'Units', 'centimeters');
    
    % change the size of the figure
    figpos = get(hfig, 'Position');
    
    set(hfig, 'Position', [figpos(1), figpos(2), figwidth, figheight]);
    
    % change the size of the paper
    set(hfig, 'PaperUnits','centimeters');
    set(hfig, 'PaperSize', [figwidth, figheight]);
    set(hfig, 'PaperPositionMode', 'manual');
    set(hfig, 'PaperPosition',[0 0 figwidth figheight]);    
    
    % reset to original units for axes and figure 
    if ~iscell(origaxunits)
        origaxunits = {origaxunits};
    end

    for i = 1:numel(hax)
        set(hax(i), 'Units', origaxunits{i});
    end

    set(hfig, 'Units', origfigunits);

function info = get_info(dname)
    cd(dname)
    files = dir('*.dcm');
    l = length(files);
    for k=1:round(size(files,1))
        info{k} = dicominfo(files(k).name);
        SliceLocation(k) = info{k}.SliceLocation;
    end
    
function [x y z x0 y0 z0 x_length y_length z_length pix_space Im_axial] = prepare_volume(V,info)
    pix_space = info{1}.PixelSpacing;
    xyz0 = info{1}.ImagePositionPatient;
    x0 = xyz0(1); y0 = xyz0(2); z0 = xyz0(3);
    x_length = size(V,2); y_length = size(V,3); z_length = size(V,1);
    x = x0 + [0:x_length-1]*pix_space(1);
    y = y0 + [0:y_length-1]*pix_space(2); 
    z = z0 - [0:z_length-1]*pix_space(1);
    for k=1:z_length
        Im_axial{k} = mat2gray(squeeze(V(k,:,:)));
        Im_axial{k} = Im_axial{k}';
    end

function [centers_desc radii_desc centers_asc radii_asc] = find_desc_asc_from_ref_scans(im,centers_desc_ref,radii_desc_ref,centers_asc_ref,radii_asc_ref);
global ObjectPolarity_mode
    I = im;
    [centers, radii] = imfindcircles(I,[1 100],'ObjectPolarity',ObjectPolarity_mode);
    %- Delete too large radii
    ind = find(radii>15);
    centers(ind,:) = [];
    radii(ind) = [];

    %% Descending aorta
    if isempty(radii_desc_ref)==1
        %- Y limits of ROI
        ymin_desc = round(size(I,1)*0.5); ymax_desc = round(size(I,1)*1);
        %- X limits of ROI
        xmin_desc = round(size(I,2)*0.4); xmax_desc = round(size(I,2)*0.6);
        %- Find circles within ROI
        centers_desc_potential = centers;radii_desc_potential = radii;
        for k=size(centers_desc_potential,1):-1:1
            if centers_desc_potential(k,1)<xmin_desc || centers_desc_potential(k,1)>xmax_desc || centers_desc_potential(k,2)<ymin_desc || centers_desc_potential(k,2)>ymax_desc
                centers_desc_potential(k,:) = [];
                radii_desc_potential(k,:) = [];
            end
        end
%     imshow(I),rectangle('Position',[xmin_desc ymin_desc xmax_desc-xmin_desc ymax_desc-ymin_desc]),viscircles(centers_desc_potential,radii_desc_potential,'Color','b')
        %- Select descending aorta
        [radii_desc imax] = max(radii_desc_potential);
        centers_desc = centers_desc_potential(imax,:);
    else
        diff_centers = []; diff_radii = []; diff_total = [];  
        for k=1:length(radii)
            diff_centers(k) = sqrt((centers(k,1)-centers_desc_ref(1))^2+(centers(k,2)-centers_desc_ref(2))^2);
            diff_radii(k) = radii(k) - radii_desc_ref;
            diff_total(k) = sqrt(diff_centers(k)^2 + diff_radii(k)^2);
        end
        [c imin] = min(diff_total);
        radii_desc = radii(imin);
        centers_desc = centers(imin,:)
    end
    
    %% Ascending aorta
    if isempty(radii_asc_ref)==1
        %- Y limits of ROI
        ymin_asc = 1; ymax_asc = round(size(I,1)*0.4);
        %- X limits of ROI
        xmin_asc = round(size(I,2)*0.25); xmax_asc = round(size(I,2)*0.5);
        %- Find circles within ROI
        centers_asc_potential = centers;radii_asc_potential = radii;
        for k=size(centers_asc_potential,1):-1:1
            if centers_asc_potential(k,1)<xmin_asc || centers_asc_potential(k,1)>xmax_asc || centers_asc_potential(k,2)<ymin_asc || centers_asc_potential(k,2)>ymax_asc
                centers_asc_potential(k,:) = [];
                radii_asc_potential(k,:) = [];
            end
        end
%     imshow(I),rectangle('Position',[xmin_asc ymin_asc xmax_asc-xmin_asc ymax_asc-ymin_asc]),viscircles(centers_asc_potential,radii_asc_potential,'Color','b')
        %- Select descending aorta
        [radii_asc imax] = max(radii_asc_potential);
        centers_asc = centers_asc_potential(imax,:);
    else
        diff_centers = []; diff_radii = []; diff_total = [];
        for k=1:length(radii)
            diff_centers(k) = sqrt((centers(k,1)-centers_asc_ref(1))^2+(centers(k,2)-centers_asc_ref(2))^2);
            diff_radii(k) = radii(k) - radii_asc_ref;
            diff_total(k) = diff_centers(k) + diff_radii(k);
        end
        [c imin] = min(diff_total);
        radii_asc = radii(imin);
        centers_asc = centers(imin,:)
    end
    
function ObjectPolarity_mode = decide_ObjectPolarity_mode(im);
    middle = [round(size(im,1)/2) round(size(im,2)/2)];
    span = [round(size(im,1)*.2) round(size(im,2)*.2)];
    im_pol = im(middle(1)-span(1):middle(1)+span(1),middle(2)-span(2):middle(2)+span(2));
    m_pol = mean(mean(im_pol));
    if m_pol >0.25
        ObjectPolarity_mode = 'bright';
    else
        ObjectPolarity_mode = 'dark';
    end
    
function [xo yo R] = circle_fit(x,y)
% A function to find the best circle fit (radius and center location) to
% given x,y pairs
% 
% Val Schmidt
% Center for Coastal and Ocean Mapping
% University of New Hampshire
% 2012
% 
% Arguments:
% x:         x coordinates
% y:         y coordinates
%
% Output:
% xo:        circle x coordinate center
% yo:        circle y coordinate center
% R:         circle radius
x = x(:);
y = y(:);
% Fitting a circle to the data - least squares style. 
%Here we start with
% (x-xo).^2 + (y-yo).^2 = R.^2
% Rewrite:
% x.^2 -2 x xo + xo^2 + y.^2 -2 y yo + yo.^2 = R.^2
% Put in matrix form:
% [-2x -2y 1 ] [xo yo -R^2+xo^2+yo^2]' = -(x.^2 + y.^2)
% Solve in least squares way...
A = [-2*x -2*y ones(length(x),1)];
x = A\-(x.^2+y.^2);
xo=x(1);
yo=x(2);
R = sqrt(  xo.^2 + yo.^2  - x(3));

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

function row = find_row(l)
    if l<8
        row =1;
    elseif l>7 && l<15
        row = 2;
    elseif l>14 && l<22
        row = 3;
    elseif l>21 && l<29
        row = 4;
    elseif l>28 && l<36
        row = 5;
    elseif l>35 && l<43
        row = 6;
    elseif l>42 && l<50
        row = 7;
    elseif l>49 && l<57
        row = 8;
    elseif l>56 && l<64
        row = 9;
    elseif l>63 && l<71
        row = 10;
    elseif l>70 && l<78
        row = 11;
    elseif l>77 && l<85
        row = 12;
    elseif l>84 && l<92
        row = 13;
    elseif l>91 && l<99
        row = 14;
    elseif l>98 && l<106
        row = 15;
    elseif l>105 && l<113
        row = 16;
    elseif l>112 && l<120
        row = 17;
    elseif l>119 && l<127
        row = 18;
    elseif l>126 && l<134
        row = 19;
    elseif l>133 && l<141
        row = 20;
    elseif l>140 && l<148
        row = 21;
    elseif l>147 && l<155
        row = 22;
    elseif l>154 && l<162
        row = 23;
    elseif l>161 && l<169
        row = 24;
    elseif l>168 && l<176
        row = 25;
    elseif l>175 && l<183
        row = 26;
    elseif l>182 && l<190
        row = 27;
    elseif l>189 && l<197
        row = 28;  
    elseif l>196 && l<204
        row = 29;
    else
        row = 30;
    end
