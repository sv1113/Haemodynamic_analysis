function varargout = GUI_Transit_Time(varargin)
% GUI_TRANSIT_TIME MATLAB code for GUI_Transit_Time.fig
%      GUI_TRANSIT_TIME, by itself, creates a new GUI_TRANSIT_TIME or raises the existing
%      singleton*.
%
%      H = GUI_TRANSIT_TIME returns the handle to a new GUI_TRANSIT_TIME or the handle to
%      the existing singleton*.
%
%      GUI_TRANSIT_TIME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_TRANSIT_TIME.M with the given input arguments.
%
%      GUI_TRANSIT_TIME('Property','Value',...) creates a new GUI_TRANSIT_TIME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Transit_Time_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Transit_Time_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Transit_Time

% Last Modified by GUIDE v2.5 09-Nov-2020 16:24:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Transit_Time_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Transit_Time_OutputFcn, ...
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


% --- Executes just before GUI_Transit_Time is made visible.
function GUI_Transit_Time_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for GUI_Transit_Time
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Transit_Time wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Transit_Time_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;
close

% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global transit_time transit_time_asc_desc
flow = getappdata(0,'flow')
freq = getappdata(0,'freq')
waveform = 2; %- velocity/flow waves
continuous = 0; %- Single wave
show = 1; %- Display wave in handles.axes1
if max(flow(:,2))<abs(min(flow(:,2))) %- Check descending flow wave is positive
    flow(:,2) = -flow(:,2);
end
switch(get(eventdata.NewValue,'Tag'));
    case 'F2F_button'
        cla, set(handles.Transit_time_str, 'String', []);
        algorithm = 1;
        transit_time = TTAlgorithm(flow,freq,algorithm,waveform,continuous,show,handles);
        TT_str = sprintf('%.5f', transit_time);
        set(handles.Transit_time_str, 'String', TT_str);
    case 'F2F_radius_button'
        cla, set(handles.Transit_time_str, 'String', []);
        algorithm = 2;
        transit_time = TTAlgorithm(flow,freq,algorithm,waveform,continuous,show,handles);
        TT_str = sprintf('%.5f', transit_time);
        set(handles.Transit_time_str, 'String', TT_str);
    case 'Least_Squares_button'    
        cla, set(handles.Transit_time_str, 'String', []);
        algorithm = 3;
        transit_time = TTAlgorithm(flow,freq,algorithm,waveform,continuous,show,handles);
        TT_str = sprintf('%.5f', transit_time);
        set(handles.Transit_time_str, 'String', TT_str);
    case 'Cross_correlation_button'
        cla, set(handles.Transit_time_str, 'String', []);
        algorithm = 4;
        transit_time = TTAlgorithm(flow,freq,algorithm,waveform,continuous,show,handles);
        TT_str = sprintf('%.5f', transit_time);
        set(handles.Transit_time_str, 'String', TT_str);
end
transit_time_asc_desc = transit_time;
%% Update flow analysis axis
Q = flow(:,1); t = [1:length(Q)]./freq;
%- Realign signal annd get timings
% [c imax] = max(Q)
% Qes = Q(1:imax)
% dQes = diff(Qes)
% [ind_intersection] = intersection(dQes)
% if isempty(ind_intersection)==0
%     Q = [Q(ind_intersection(end):end) ; Q(1:ind_intersection(end)-1)]
% end
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
setappdata(0,'EF1',EF1)

%- HR
HR = 60/t(end);
HR_str = sprintf('%.1f', HR);
set(handles.str_HR, 'String', HR_str);
setappdata(0,'HR',HR)
%- SV
SV = trapz(t,Q);
SV_str = sprintf('%.1f', SV);
set(handles.str_SV, 'String', SV_str);
setappdata(0,'SV',SV)
%- CO
CO = SV*1e-3*HR;
CO_str = sprintf('%.1f', CO);
set(handles.str_CO, 'String', CO_str);
setappdata(0,'CO',CO)
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in Save_button.
function Save_button_Callback(hObject, eventdata, handles)
global transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia length_asc_desc length_desc_dia length_asc_dia ...
        transit_time path_length PWV Tasc Tdesc Tdia Qasc Qdesc Qdia Aasc Adesc Adia PWV_asc_desc PWV_desc_dia PWV_asc_dia ...
        info Dis_Asc Dis_Desc Dis_Dia Tasc Tdesc Tdia Qasc Qdesc Qdia Aasc Adesc Adia
setappdata(0,'transit_time_asc_desc',transit_time_asc_desc)
setappdata(0,'transit_time_desc_dia',transit_time_desc_dia)
setappdata(0,'transit_time_asc_dia',transit_time_asc_dia)
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
    [file,path,indx] = uiputfile('.xls');
    generate_spreadsheet_1_part(file,path,info,EF1,SV,HR,CO,transit_time_asc_desc,PWV_asc_desc,Dis_Asc,Dis_Desc,Tasc,Qasc,Qdesc,Aasc,Adesc)
else    
    Time_s = Tasc; Flow_Ascending_ml_s = Qasc; Flow_Descending_ml_s = Qdesc; Flow_Diaphragm_ml_s = Qdia;
    Area_Ascending_mm2 = Aasc; Area_Descending_mm2 = Adesc; Area_Diaphragm_mm2 = Adia; 
    TT_Ascending_Descending_s = transit_time_asc_desc; TT_Descending_Diaphragm_s = transit_time_desc_dia; TT_Ascending_Diaphragm_s = transit_time_asc_dia;
    PWV_Ascending_Descending_m_s = PWV_asc_desc; PWV_Descending_Diaphragm_m_s = PWV_desc_dia; PWV_Ascending_Diaphragm_m_s = PWV_asc_dia;
    EF1_pct = EF1; SV_ml = SV; HR_bpm = HR; CO_m3_s = CO;
    cd(path_flow);
    [file,path,indx] = uiputfile('.xls');
    generate_spreadsheet_3_parts(file,path,info,EF1,SV,HR,CO,transit_time_asc_desc,transit_time_desc_dia,transit_time_asc_dia,PWV_asc_desc,PWV_desc_dia,PWV_asc_dia,Dis_Asc,Dis_Desc,Dis_Dia,Tasc,Qasc,Qdesc,Qdia,Aasc,Adesc,Adia)
end
cd(path_origin)
uiresume

% --- Executes on button press in Continue_button.
function Continue_button_Callback(hObject, eventdata, handles)
% hObject    handle to Continue_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global transit_time
setappdata(0,'transit_time_asc_desc',transit_time)
uiresume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADDITIONAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TT = TTAlgorithm(signal,f,algorithm,waveform,continuous,show,handles)
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
[indmaxs,indmins,gradins] = Find_Landmarks(t,signal,t_intv,npoly_1,npoly1_1,f,algorithm,waveform,continuous,show,handles);

switch algorithm
    case 3
        % Least Squares difference of the systolic upstroke
        t_int = PWV_Calulator_Least_Squares(t,signal,indmaxs,indmins,handles);
    case 4
        % Cross correlation of the entire waveform
        t_int = PWV_Calulator_CC_cycle(t,signal,indmaxs,indmins,gradins,show,handles);
    otherwise
        % Both foot to foot algorithms
        t_int = PWV_Calulator_FTF(t,signal,indmaxs,indmins,gradins,show,handles);
end
TT = t_int;

function [t_int,t_foot,horint] = PWV_Calulator_FTF(t,signal,indmaxs,indmins,gradins,show,handles)

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
        axes(handles.axes1)
        plot(t_foot,horint,'rx','MarkerSize',12,'Linewidth',2)
end

non_zeros = find(t_foot(1,:).*t_foot(2,:));

t_int = TT;

function t_int= PWV_Calulator_Least_Squares(t,signal,indmaxs,indmins,handles)

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

function t_int = PWV_Calulator_CC_cycle(t,signal,indmaxs,indmins,gradins,show,handles)
signal = signal;
t = t;

dt = t(2)-t(1);

shift = round(0.2/dt);

TTApprox = median(diff(gradins(:,:,1)));

axes(handles.axes1)
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

function [indmaxs,indmins,gradins] = Find_Landmarks(t,signal,t_intv,npoly_1,npoly1_1,f,algorithm,waveform,continuous,show,handles)
%%%%%%%%%%%%% Locate Minima, Maxima and Max Gradients %%%%%%%%%%%%%
t = t';
signal = signal';

switch show
    case 0
    case 1
        axes(handles.axes1)
        cla
        hold on
        p1 = plot(t,signal(:,1),'b','linewidth',1);
        p2 = plot(t,signal(:,2),'r','linewidth',1);
        legend([p1 p2],{'Ascending','Descending'},'AutoUpdate','off')
        legend boxoff
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
            axes(handles.axes1)
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

%%%%%%%%%%%%% Determinsetappdata(0,'time_transit',transit_time)e Minima, Maxima and Max Gradients - End %%%%%%%%%%%%%%%%%


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
        axes(handles.axes1)
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
        test = flow(k-1)*flow(k)
        if test<0
            ind_intersection(m) = k;
            m=m+1;
        end
    end

function [] = generate_spreadsheet_1_part(filename,path,info,EF1,SV,HR,CO,transit_time_asc_desc,PWV_asc_desc,Dis_Asc,Dis_Desc,Tasc,Qasc,Qdesc,Aasc,Adesc)

cd(path)

%% Global info
filename = [filename,'.xls'];
sheet=1;

%% Patient Data
xlRange='A1';
Results_Names={'%%% PATIENT INFO'};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A2';
Results_Names={'Name',[info.PatientName.GivenName,' ',info.PatientName.FamilyName]};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A3';
Results_Names={'ID',[info.PatientID]};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A4';
Results_Names={'Gender',[info.PatientSex]};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A5';
DOB = [num2str(info.PatientBirthDate(7:8)),'/',num2str(info.PatientBirthDate(5:6)),'/',num2str(info.PatientBirthDate(1:4))];
Results_Names={'Date Of Birth',DOB};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A6';
Date_Scan = [num2str(info.StudyDate(7:8)),'/',num2str(info.StudyDate(5:6)),'/',num2str(info.StudyDate(1:4))];
Results_Names={'Date Of Scan',Date_Scan};
xlswrite(filename,Results_Names,sheet,xlRange);

%% Flow Analysis
xlRange='A8';
Results_Names={'Flow analysis'};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A9';
Results_Names={'EF1 (%):',[num2str(EF1,'%.2f\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A10';
Results_Names={'SV (ml):',[num2str(SV,'%.2f\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A11';
Results_Names={'HR (bpm):',[num2str(HR,'%.2f\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A12';
Results_Names={'CO (L/min):',[num2str(CO,'%.2f\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);

%% PWV analysis
%- Ascending-Descending
xlRange='A14';
Results_Names={'Ascending-Descending PWV'};
xlswrite(filename,Results_Names,sheet,xlRange);
% xlRange='A15';
% Results_Names={'Path length (m):',[num2str(length_asc_desc/1000,'%.3f\n')]};
% xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A16';
Results_Names={'Transit time (s):',[num2str(transit_time_asc_desc,'%.3f\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);
% xlRange='A17';
% Results_Names={'PWV (m/s):',[num2str(PWV_asc_desc,'%.2f\n')]};
% xlswrite(filename,Results_Names,sheet,xlRange);

%% Distensibility analysis
xlRange='A18';
Results_Names={'%%% DISTENSIBILITY'};
xlswrite(filename,Results_Names,sheet,xlRange);
% xlRange='A19';
% Results_Names={'Ascending Aorta Distensibility (/mmHg):',[num2str(Dis_Asc,'%.2e\n')]};
% xlswrite(filename,Results_Names,sheet,xlRange);

%% Plot Time, Flow and Area waveforms
xlRange='A21';
Results_Names={'Waveforms'};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A22';
Results_Names={'Time (ms)','Q ascending (ml/s)','Q Descending (ml/s)'};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A23';
Results_Values=[Tasc Qasc Qdesc];
xlswrite(filename,Results_Values,sheet,xlRange);

function [] = generate_spreadsheet_3_parts(filename,path,info,EF1,SV,HR,CO,transit_time_asc_desc,transit_time_desc_dia,transit_time_asc_dia,PWV_asc_desc,PWV_desc_dia,PWV_asc_dia,Dis_Asc,Dis_Desc,Dis_Dia,Tasc,Qasc,Qdesc,Qdia,Aasc,Adesc,Adia)

cd(path)

%% Global info
sheet=1;

%% Patient Data
xlRange='A1';
Results_Names={'%%% PATIENT INFO'};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A2';
Results_Names={'Name',[info.PatientName.GivenName,' ',info.PatientName.FamilyName]};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A3';
Results_Names={'ID',[info.PatientID]};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A4';
Results_Names={'Gender',[info.PatientSex]};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A5';
DOB = [num2str(info.PatientBirthDate(7:8)),'/',num2str(info.PatientBirthDate(5:6)),'/',num2str(info.PatientBirthDate(1:4))];
Results_Names={'Date Of Birth',DOB};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A6';
Date_Scan = [num2str(info.StudyDate(7:8)),'/',num2str(info.StudyDate(5:6)),'/',num2str(info.StudyDate(1:4))];
Results_Names={'Date Of Scan',Date_Scan};
xlswrite(filename,Results_Names,sheet,xlRange);

%% Flow Analysis
xlRange='A8';
Results_Names={'%%% FLOW ANALYSIS'};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A9';
Results_Names={'EF1 (%):',[num2str(EF1,'%.2f\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A10';
Results_Names={'SV (ml):',[num2str(SV,'%.2f\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A11';
Results_Names={'HR (bpm):',[num2str(HR,'%.2f\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A12';
Results_Names={'CO (L/min):',[num2str(CO,'%.2f\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);

%% PWV analysis
%- Ascending-Descending
xlRange='A14';
Results_Names={'%%% ASCENDING-DESCENDING PWV'};
xlswrite(filename,Results_Names,sheet,xlRange);
% xlRange='A15';
% Results_Names={'Path length (m):',[num2str(length_asc_desc/1000,'%.3f\n')]};
% xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A16';
Results_Names={'Transit time (s):',[num2str(transit_time_asc_desc,'%.3f\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);
% xlRange='A17';
% Results_Names={'PWV (m/s):',[num2str(PWV_asc_desc,'%.2f\n')]};
% xlswrite(filename,Results_Names,sheet,xlRange);

%- Descending-Diaphragm
xlRange='A19';
Results_Names={'%%% DESCENDINIG-DIAPHRAGM PWV'};
xlswrite(filename,Results_Names,sheet,xlRange);
% xlRange='A20';
% Results_Names={'Path length (m):',[num2str(length_desc_dia/1000,'%.3f\n')]};
% xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A21';
Results_Names={'Transit time (s):',[num2str(transit_time_desc_dia,'%.3f\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);
% xlRange='A22';
% Results_Names={'PWV (m/s):',[num2str(PWV_desc_dia,'%.2f\n')]};
% xlswrite(filename,Results_Names,sheet,xlRange);

%- Ascending-Diapragm
xlRange='A24';
Results_Names={'%%% ASCENDING-DIAPHRAGM PWV'};
xlswrite(filename,Results_Names,sheet,xlRange);
% xlRange='A25';
% Results_Names={'Path length (m):',[num2str(length_asc_dia/1000,'%.3f\n')]};
% xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A26';
Results_Names={'Transit time (s):',[num2str(transit_time_asc_dia,'%.3f\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);
% xlRange='A27';
% Results_Names={'PWV (m/s):',[num2str(PWV_asc_dia,'%.2f\n')]};
% xlswrite(filename,Results_Names,sheet,xlRange);

%% Distensibility analysis
xlRange='A29';
Results_Names={'%%% DISTENSIBILITY'};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A30';
Results_Names={'Ascending Aorta Distensibility (/mmHg):',[num2str(Dis_Asc,'%.2e\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A31';
Results_Names={'Descending Aorta Distensibility (/mmHg):',[num2str(Dis_Desc,'%.2e\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);
xlRange='A32';
Results_Names={'Diaphragm Aorta Distensibility (/mmHg):',[num2str(Dis_Dia,'%.2e\n')]};
xlswrite(filename,Results_Names,sheet,xlRange);

%% Plot Time, Flow and Area waveforms
xlRange='A34';
Results_Names={'%%% WAVEFORMS'};
xlswrite(filename,Results_Names,sheet,xlRange);
if ~isempty(Aasc)
    xlRange='A35';
    Results_Names={'Time (s)','Q ascending (ml/s)','Q descending (ml/s)','Q diaphragm (ml/s)','A ascending (mm2)','A descending (mm2)','A diaphragm (mm2)'};
    xlswrite(filename,Results_Names,sheet,xlRange);
    xlRange='A36';
    Results_Values=[Tasc Qasc Qdesc Qdia Aasc Adesc Adia];
    xlswrite(filename,Results_Values,sheet,xlRange);
else
    xlRange='A35';
    Results_Names={'Time (ms)','Q ascending (ml/s)','Q descending (ml/s)','Q diaphragm (ml/s)'};
    xlswrite(filename,Results_Names,sheet,xlRange);
    xlRange='A36';
    Results_Values=[Tasc Qasc Qdesc Qdia];
    xlswrite(filename,Results_Values,sheet,xlRange);
end

    