function varargout = GUI_PWV_coronal_manual_choice(varargin)
% GUI_PWV_CORONAL_MANUAL_CHOICE MATLAB code for GUI_PWV_manual_coronal.fig
%      GUI_PWV_CORONAL_MANUAL_CHOICE, by itself, creates a new GUI_PWV_CORONAL_MANUAL_CHOICE or raises the existing
%      singleton*.
%
%      H = GUI_PWV_CORONAL_MANUAL_CHOICE returns the handle to a new GUI_PWV_CORONAL_MANUAL_CHOICE or the handle to
%      the existing singleton*.
%
%      GUI_PWV_CORONAL_MANUAL_CHOICE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PWV_CORONAL_MANUAL_CHOICE.M with the given input arguments.
%
%      GUI_PWV_CORONAL_MANUAL_CHOICE('Property','Value',...) creates a new GUI_PWV_CORONAL_MANUAL_CHOICE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_PWV_manual_auto_choiceOpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_PWV_manual_auto_choice_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help
% GUI_PWV_manual_auto_choice

% Last Modified by GUIDE v2.5 21-Apr-2021 15:44:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_PWV_coronal_manual_choice_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_PWV_coronal_manual_choice_OutputFcn, ...
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


% --- Executes just before GUI_PWV_manual_auto_choice is made visible.
function GUI_PWV_coronal_manual_choice_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_PWV_coronal_manual_choice (see VARARGIN)

% Choose default command line output for GUI_PWV_coronal_manual_choice
handles.output = hObject;

% Initiate and clear variables
global transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia info info_global C_total...
    length_asc_desc length_desc_dia length_asc_dia PWV_asc_desc PWV_desc_dia PWV_asc_dia row col...
    ObjectPolarity_mode flag_diapgragm z_ref z_ref_diaphragm
axes(handles.axes1); cla reset;

%- Retrieve coordinates and radius of reference ascending/descending slice
outputname = {'centers_desc','radii_desc','centers_asc','radii_asc'};
k_ref_asc = getappdata(0,'k_ref_asc')
for k=1:length(outputname)
    eval([outputname{k} '(k_ref_asc,:)=getappdata(0,''' outputname{k} ''');'])
end
centers_desc_ref = centers_desc(k_ref_asc,:),radii_desc_ref = radii_desc(k_ref_asc),centers_asc_ref = centers_asc(k_ref_asc,:),radii_asc_ref = radii_asc(k_ref_asc,:),
%- Retrieve parameters
outputname = {'transit_time_asc_desc','transit_time_desc_dia','transit_time_asc_dia','ObjectPolarity_mode','row','col',...
    'info','info_global','Im_axial','x','y','z','k_ref_asc_plot','k_ref_dia','V','flag_diapgragm','pix_space','Dis_Asc','Dis_Desc',...
    'Dis_Dia','EDV','EF','DBP','SBP','Pes','ENd_est','Ees','PP','z_ref','z_ref_diaphragm'};
for k=1:length(outputname)
    eval([outputname{k} '=getappdata(0,''' outputname{k} ''');'])
end

setappdata(0,'length_asc_desc',[]),setappdata(0,'length_desc_dia',[]),setappdata(0,'length_asc_dia',[]),
setappdata(0,'PWV_asc_desc',[]),setappdata(0,'PWV_desc_dia',[]),setappdata(0,'PWV_asc_dia',[]),
set(handles.str_PWV_asc_desc,'String','');set(handles.str_PWV_desc_dia,'String','');set(handles.str_PWV_asc_dia,'String','');
set(handles.str_length_asc_desc,'String','');set(handles.str_length_desc_dia,'String','');set(handles.str_length_asc_dia,'String','');
set(handles.str_TT_asc_desc,'String',round(transit_time_asc_desc,3));set(handles.str_TT_desc_dia,'String',round(transit_time_desc_dia,3));set(handles.str_TT_asc_dia,'String',round(transit_time_asc_dia,3));

length_asc_desc = []; length_desc_dia = []; length_asc_dia = [];
PWV_asc_desc = []; PWV_desc_dia = []; PWV_asc_dia = [];
% Update handles structure
guidata(hObject, handles);

%- Plot ascending and descending aorta for reference slice
k_plot = k_ref_asc_plot;
for k=1
    teta=0:0.01:2*pi ; flag_plot = 0;
    axes(handles.axes1),hold all,view([107 17]),
%     figure(1),hold all,
    %- Display reference plane
    x_min = min(centers_desc_ref(1)-radii_desc_ref,centers_asc_ref(1)-radii_asc_ref)-10;x_max = max(centers_desc_ref(1)+radii_desc_ref,centers_asc_ref(1)+radii_asc_ref)+10;
    y_min = min(centers_desc_ref(2)-radii_desc_ref,centers_asc_ref(2)-radii_asc_ref)-10;y_max = max(centers_desc_ref(2)+radii_desc_ref,centers_asc_ref(2)+radii_asc_ref)+10;
%     x_plane = [x_min x_max ; x_min x_max]
    y_plane = [1 1 ; size(Im_axial{k_ref_asc},1) size(Im_axial{k_ref_asc},1)]
    x_plane = [1 size(Im_axial{k_ref_asc},2) ; 1 size(Im_axial{k_ref_asc},2)]
%     y_plane = [y_min y_min ; y_max y_max]
    z_plane_asc_desc = [z(k_ref_asc) z(k_ref_asc) ; z(k_ref_asc) z(k_ref_asc)]
    axes(handles.axes1), 
    surf(x_plane,y_plane,z_plane_asc_desc,...    % Plot the surface
     'CData',Im_axial{k_ref_asc},...
     'FaceAlpha', 0.6,...
     'FaceColor','texturemap'); 
%      'FaceColor','texturemap');
    %- Display slice where diaphragm flow was acquired
    if flag_diapgragm == 1
        z_plane_diaphragm = [z(k_ref_dia) z(k_ref_dia) ; z(k_ref_dia) z(k_ref_dia)]
        surf(x_plane,y_plane,z_plane_diaphragm,...    % Plot the surface
        'CData',Im_axial{k_ref_dia},...
        'FaceAlpha', 0.6,...
        'FaceColor','texturemap');
        colormap(gray)
    end
   %- 3D plot
    %- Descending aorta
    C_desc_3D(k_ref_asc,:) = [centers_desc(k_ref_asc,:) z(k_ref_asc)];
    R_desc_3D(k_ref_asc) = radii_desc(k_ref_asc);
    x_desc = C_desc_3D(k_ref_asc,1)+R_desc_3D(k_ref_asc)*cos(teta);
    y_desc = C_desc_3D(k_ref_asc,2)+R_desc_3D(k_ref_asc)*sin(teta) ;
    z_desc = C_desc_3D(k_ref_asc,3)+zeros(size(x_desc)) ;
    axes(handles.axes1),patch(x_desc,y_desc,z_desc,'b');
    axes(handles.axes1),
%     plot3(C_desc_3D(k_ref_asc,1),C_desc_3D(k_ref_asc,2),C_desc_3D(k_ref_asc,3),'*k'); 
    %- Ascending aorta
    C_asc_3D(k_ref_asc,:) = [centers_asc(k_ref_asc,:) z(k_ref_asc)];
    R_asc_3D(k_ref_asc) = radii_asc(k_ref_asc);
    x_asc = C_asc_3D(k_ref_asc,1)+R_asc_3D(k_ref_asc)*cos(teta);
    y_asc = C_asc_3D(k_ref_asc,2)+R_asc_3D(k_ref_asc)*sin(teta) ;
    z_asc = C_asc_3D(k_ref_asc,3)+zeros(size(x_desc)) ;
    axes(handles.axes1),patch(x_asc,y_asc,z_asc,'r');
    axes(handles.axes1),
%     plot3(C_asc_3D(k_ref_asc,1),C_asc_3D(k_ref_asc,2),C_asc_3D(k_ref_asc,3),'*k'); 
end

%- Find ascending and descending aorta for slices below reference slice
flag_asc = 1; 
if flag_diapgragm==1
    for k=k_ref_asc+1:5:k_ref_dia+5
        Im = Im_axial{k};
        setappdata(0,'im',Im)
        
        %- Set up ascending and descending aorta from previous slices
        setappdata(0,'centers_desc_ref',centers_desc_ref)
        setappdata(0,'radii_desc_ref',radii_desc_ref)
        setappdata(0,'centers_asc_ref',centers_asc_ref)
        setappdata(0,'radii_asc_ref',radii_asc_ref)
        GUI_ref_asc_desc_manual
        
        %- Retrieve information about reference sections from first GUI
        outputname = {'centers_desc','radii_desc','centers_asc','radii_asc'};
        for j=1:length(outputname)
            eval([outputname{j} '(k,:)=getappdata(0,''' outputname{j} ''');'])
        end
        centers_desc_ref = centers_desc(k,:); radii_desc_ref = radii_desc(k);
        centers_asc_ref = centers_asc(k,:); radii_asc_ref = radii_asc(k);
        
        %- Plot the aortic sections
        C_desc_3D(k,:) = [centers_desc(k,:) z(k)]; R_desc_3D(k) = radii_desc(k);
        C_asc_3D(k,:) = [centers_asc(k,:) z(k)]; R_asc_3D(k) = radii_asc(k);
        %- Define descending aorta for 3D plot
        x_desc = C_desc_3D(k,1)+R_desc_3D(k)*cos(teta); y_desc = C_desc_3D(k,2)+R_desc_3D(k)*sin(teta) ;
        z_desc = C_desc_3D(k,3)+zeros(size(x_desc));
        axes(handles.axes1),patch(x_desc,y_desc,z_desc,'b'); axes(handles.axes1),
%         plot3(C_desc_3D(k,1),C_desc_3D(k,2),C_desc_3D(k,3),'*k'); 
        %- Define ascending aorta for 3D plot
        if radii_asc(k)~=0
            x_asc = C_asc_3D(k,1)+R_asc_3D(k)*cos(teta); y_asc = C_asc_3D(k,2)+R_asc_3D(k)*sin(teta) ;
            z_asc = C_asc_3D(k,3)+zeros(size(x_asc));
            axes(handles.axes1),patch(x_asc,y_asc,z_asc,'r'); axes(handles.axes1),
%             plot3(C_asc_3D(k,1),C_asc_3D(k,2),C_asc_3D(k,3),'*k');
            flag_asc = 1;
        else
            flag_asc = 0;
        end
        setappdata(0,'flag_asc',flag_asc)
    end
end


%- Find ascending and descending aorta for slices above reference slice
setappdata(0,'flag_asc',1)
k_plot = k_ref_asc_plot+1;
centers_desc_ref = centers_desc(k_ref_asc,:); radii_desc_ref = radii_desc(k_ref_asc); centers_asc_ref = centers_asc(k_ref_asc,:); radii_asc_ref = radii_asc(k_ref_asc);
for k= k_ref_asc-1:-3:k_ref_asc-100
        Im = Im_axial{k};
        setappdata(0,'im',Im)
        
        %- Set up ascending and descending aorta from previous slices
        setappdata(0,'centers_desc_ref',centers_desc_ref)
        setappdata(0,'radii_desc_ref',radii_desc_ref)
        setappdata(0,'centers_asc_ref',centers_asc_ref)
        setappdata(0,'radii_asc_ref',radii_asc_ref)
        GUI_ref_asc_desc_manual
        
        %- Retrieve information about reference sections from first GUI
        flag_quit = getappdata(0,'flag_quit')
        if flag_quit==0
            outputname = {'centers_desc','radii_desc','centers_asc','radii_asc'};
            for j=1:length(outputname)
                eval([outputname{j} '(k,:)=getappdata(0,''' outputname{j} ''');'])
            end
            centers_desc_ref = centers_desc(k,:); radii_desc_ref = radii_desc(k);
            centers_asc_ref = centers_asc(k,:); radii_asc_ref = radii_asc(k);
        
            C_desc_3D(k,:) = [centers_desc(k,:) z(k)]; R_desc_3D(k) = radii_desc(k);
            C_asc_3D(k,:) = [centers_asc(k,:) z(k)]; R_asc_3D(k) = radii_asc(k);
            %- Define descending aorta for 3D plot
            x_desc = C_desc_3D(k,1)+R_desc_3D(k)*cos(teta); y_desc = C_desc_3D(k,2)+R_desc_3D(k)*sin(teta) ;
            z_desc = C_desc_3D(k,3)+zeros(size(x_desc));
            axes(handles.axes1),patch(x_desc,y_desc,z_desc,'b'); axes(handles.axes1),
%             plot3(C_desc_3D(k,1),C_desc_3D(k,2),C_desc_3D(k,3),'*k'); 
            %- Define ascending aorta for 3D plot
            if radii_asc(k)~=0
                x_asc = C_asc_3D(k,1)+R_asc_3D(k)*cos(teta); y_asc = C_asc_3D(k,2)+R_asc_3D(k)*sin(teta) ;
                z_asc = C_asc_3D(k,3)+zeros(size(x_asc));
                axes(handles.axes1),patch(x_asc,y_asc,z_asc,'r'); axes(handles.axes1),
%                 plot3(C_asc_3D(k,1),C_asc_3D(k,2),C_asc_3D(k,3),'*k');
                flag_asc = 1;
            else
                flag_asc = 0;
            end
            setappdata(0,'flag_asc',flag_asc)
        else
            break
        end
end

%- Find top of the aorta
k_ref_asc = k_ref_asc;
Vtop = V(1:k_ref_asc,:,:);
ztop = z(1:k_ref_asc); xtop = x; ytop = y;
%- Find middle slice (i_ref_plot)
ind_yasc = find(C_asc_3D(:,2)>0);
ind_ydesc = find(C_desc_3D(:,2)>0);
i_ref_plot = round(mean([C_asc_3D(ind_yasc(1),2) C_desc_3D(ind_ydesc(1),2)]))
Im = (squeeze(Vtop(:,:,i_ref_plot)));
Im = mat2gray(Im);
setappdata(0,'Im_top',Im)
%- Extract top of the aorta from slice
GUI_ref_top_aorta
centers_top = getappdata(0,'centers_top')
radii_top = getappdata(0,'radii_top')
C_ref_orignal = centers_top;
C_ref = centers_top;
radii_ref = radii_top;
x3D = C_ref(1);
y3D = (i_ref_plot); 
z3D = z_plane_asc_desc(1) + (size(Vtop,1)-C_ref(2))*pix_space(1);
C_3D(i_ref_plot,:) = [x3D y3D z3D]

%- Stitch parts together
ind = find(C_asc_3D(:,1)~=0);
C1x = flipud(C_asc_3D(ind,1)); C1y =  flipud(C_asc_3D(ind,2)); C1z = flipud(C_asc_3D(ind,3));
C1 = [C1x C1y C1z];
C2x = C_3D(C_3D(:,1)~=0,1); C2y = C_3D(C_3D(:,2)~=0,2); C2z = C_3D(C_3D(:,3)~=0,3);
C2 = [C2x C2y C2z];
C3x = C_desc_3D(C_desc_3D(:,1)~=0,1); C3y = C_desc_3D(C_desc_3D(:,2)~=0,2); C3z = C_desc_3D(C_desc_3D(:,3)~=0,3);
C3 = [C3x C3y C3z];
C_total = [C1 ; C2 ; C3];
Q
%- Interpolate signal
C_total = interparc(100,C_total(:,1),C_total(:,2),C_total(:,3),'spline');

%- Plot spline
axes(handles.axes1),plot3(C_total(:,1),C_total(:,2),C_total(:,3),'-*k')

%- Calculate aortic_descending path length
ind_top = find(C_total(:,3)>z_ref)
Xtop = C_total(ind_top,1); Ytop = C_total(ind_top,2); Ztop = C_total(ind_top,3); 
length_asc_desc = calculate_aortic_arch_3D(Xtop,Ytop,Ztop);
length_asc_desc = length_asc_desc*info{1}.PixelSpacing(1);
str_length_asc_desc = sprintf('%.3f', length_asc_desc/1e3);
set(handles.str_length_asc_desc, 'String', str_length_asc_desc);

%- Calculate length from descending aorta down to diaphragm
if flag_diapgragm == 1
    ind_down = find(C_total(:,3)<z_ref & C_total(:,3)>z_ref_diaphragm);
    Xdown = C_total(ind_down,1); Ydown = C_total(ind_down,2); Zdown = C_total(ind_down,3); 
    length_desc_dia = calculate_aortic_arch_3D(Xdown,Ydown,Zdown);
    length_desc_dia = length_desc_dia*info{1}.PixelSpacing(1);
    str_length_desc_dia = sprintf('%.3f', length_desc_dia/1e3);
    set(handles.str_length_desc_dia, 'String', str_length_desc_dia);
end


%- Write PWV value if transit time value is known
if ~isempty(transit_time_asc_desc)==1
    PWV_asc_desc = length_asc_desc*1e-3/transit_time_asc_desc;
    pwv_asc_desc_str = sprintf('%.1f', PWV_asc_desc);
    set(handles.str_PWV_asc_desc, 'String', pwv_asc_desc_str);
end
if ~isempty(transit_time_desc_dia)==1
    PWV_desc_dia = length_desc_dia*1e-3/transit_time_desc_dia;
    pwv_desc_dia_str = sprintf('%.1f', PWV_desc_dia);
    set(handles.str_PWV_desc_dia, 'String', pwv_desc_dia_str);
end
if ~isempty(transit_time_asc_dia)==1
    length_asc_dia = length_asc_desc + length_desc_dia
    str_length_asc_dia = sprintf('%.3f', length_asc_dia/1e3);
    set(handles.str_length_asc_dia, 'String', str_length_asc_dia);
    PWV_asc_dia = length_asc_dia*1e-3/transit_time_asc_dia;
    pwv_asc_dia_str = sprintf('%.1f', PWV_asc_dia);
    set(handles.str_PWV_asc_dia, 'String', pwv_asc_dia_str);
end

% cd(current_folder)

guidata(hObject,handles);



% UIWAIT makes GUI_PWV_coronal_manual_choice wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_PWV_coronal_manual_choice_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in Generate_spreadsheet_button.
function Generate_spreadsheet_button_Callback(hObject, eventdata, handles)
global transit_time_asc_desc transit_time_desc_dia transit_time_asc_dia length_asc_desc length_desc_dia length_asc_dia ...
        Tasc Tdesc Tdia Qasc Qdesc Qdia Aasc Adesc Adia PWV_asc_desc PWV_desc_dia PWV_asc_dia ...
        info Dis_Asc Dis_Desc Dis_Dia Tasc Tdesc Tdia Qasc Qdesc Qdia Aasc Adesc Adia path_flow info_global C_total
EF1 = getappdata(0,'EF1'); HR = getappdata(0,'HR'); SV = getappdata(0,'SV'); CO = getappdata(0,'CO');
EDV = getappdata(0,'EDV'); EF = getappdata(0,'EF'); PP = getappdata(0,'PP'); DBP = getappdata(0,'DBP'); SBP = getappdata(0,'SBP');
Pes = getappdata(0,'Pes'); ENd_est = getappdata(0,'ENd_est'); Ees = getappdata(0,'Ees'); tNd = getappdata(0,'tNd'); ENd_avg = getappdata(0,'ENd_avg');
Dis_Asc = getappdata(0,'Dis_Asc'); Dis_Desc = getappdata(0,'Dis_Desc'); Dis_Dia = getappdata(0,'Dis_Dia'); 

path_origin = cd;
cd(path_flow);
if isempty(Tdia)
    [file,path,indx] = uiputfile('.xls');;
    generate_spreadsheet_1_part(file,path,info_global,EF1,SV,HR,CO,length_asc_desc,transit_time_asc_desc,PWV_asc_desc,Dis_Asc,Dis_Desc,Tasc,Qasc,Qdesc,Aasc,Adesc,EDV,EF,PP,DBP,SBP,Pes,ENd_est,Ees,tNd,ENd_avg,C_total)
else    
    [file,path,indx] = uiputfile('.xls');;
    generate_spreadsheet_3_parts(file,path,info_global,EF1,SV,HR,CO,length_asc_desc,length_desc_dia,length_asc_dia,transit_time_asc_desc,transit_time_desc_dia,transit_time_asc_dia,PWV_asc_desc,PWV_desc_dia,PWV_asc_dia,Dis_Asc,Dis_Desc,Dis_Dia,Tasc,Qasc,Qdesc,Qdia,Aasc,Adesc,Adia,EDV,EF,PP,DBP,SBP,Pes,ENd_est,Ees,tNd,ENd_avg,C_total)
end
cd(path_origin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADDITIONAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X Y Z] = extrapolate_top_candy_cane(centers_asc,centers_desc,SliceLocation,k_ref)
    %% Prepare vectors
    ind = find(centers_asc(:,1)~=0);
    x = [centers_asc(ind,1);flip(centers_desc(:,1))];
    y = [centers_asc(ind,2);flip(centers_desc(:,2))];
    z = [SliceLocation(ind)';flip(SliceLocation(1:size(centers_desc,1)))'];
%     plot3(x,y,z,'-*r'),view([90 0])
    %- Prepare interpolation boundaries
    y_1 = centers_asc(ind(end),2); y_2 = centers_desc(end,2);
    x_1 = centers_asc(ind(end),1); x_2 = centers_desc(end,1);
    
    %% Retrieve points above reference scan
    z_ref = SliceLocation(k_ref);
    ind = find(z>=z_ref);
    x_top = x(ind);
    y_top = y(ind);
    z_top = z(ind);
%     plot3(x_top,y_top,z_top,'-*r'),view([90 0])

    %% Interpolate data
    l = (min(y_top)+max(y_top))/20;
    y_query = [min(y_top)-l:2:max(y_top)+l];
    coef = polyfit(y_top,z_top,3);
    z_poly = polyval(coef,y_query);
    y_interpolate = y_1:1:y_2;
    z_interpolate = interp1(y_query,z_poly,y_interpolate)
    lx = (x_2-x_1)/(length(y_interpolate)-1)
    x_interpolate = x_1:lx:x_2;
%     %- 2D
%     figure(),hold all,
%     plot(y_top,z_top,'*r');
%     plot(y_interpolate,z_interpolate,'*k');
%     plot(y_query,z_poly,'g');
%     %- 3D
%     figure(),plot3(x_interpolate,y_interpolate,z_interpolate)
    
    %% Prepare output vectors
    ind = find(centers_asc(:,1)~=0); ind(end:end) = [];
    X = [centers_asc(ind,1);x_interpolate';flip(centers_desc(1:end-1,1))];
    Y = [centers_asc(ind,2);y_interpolate';flip(centers_desc(1:end-1,2))];
    Z = [ SliceLocation(ind)';z_interpolate';flip(SliceLocation(1:size(centers_desc,1)-1))'];

function distance = distance_from_reference(list,ref)
    l = size(list,1);
    for k=1:l
        x_dist = list(k,1)-ref(1);
        y_dist = list(k,2)-ref(2);
        distance(k) = sqrt(x_dist^2 + y_dist^2);
    end
    
function L = calculate_aortic_arch_3D(X_path,Y_path,Z_path);
    for k=2:length(X_path)
        X_dist = X_path(k)-X_path(k-1);
        Y_dist = Y_path(k)-Y_path(k-1);
        Z_dist = Z_path(k)-Z_path(k-1);
        l_ind(k-1) = sqrt(X_dist^2 + Y_dist^2 + Z_dist^2);
    end
    L = sum(l_ind);
    
function [centers_desc_list radii_desc_list centers_asc_list radii_asc_list] = select_potential_desc_asc(I,centers,centers_desc_ref,centers_asc_ref,radii,flag_plot)
    
    %% Descending Aorta
    %- Define rectangle
    y_width = size(I,1)/12.5; x_width = size(I,2)/12.5;
    xmin_desc = centers_desc_ref(1)-x_width; xmax_desc = centers_desc_ref(1)+x_width;
    ymin_desc = centers_desc_ref(2)-y_width; ymax_desc = centers_desc_ref(2)+x_width;
    %- Find potential descending aorta within rectangle
    centers_desc_list = centers; radii_desc_list = radii;
    for k=size(centers_desc_list,1):-1:1
        if centers_desc_list(k,1)<xmin_desc || centers_desc_list(k,1)>xmax_desc || centers_desc_list(k,2)<ymin_desc || centers_desc_list(k,2)>ymax_desc
            centers_desc_list(k,:) = [];
            radii_desc_list(k) = [];
        end
    end
    %- Delete too large circles
    ind = find(radii_desc_list>15)
    radii_desc_list(ind) = [];
    centers_desc_list(ind,:) = [];
%     figure(),imshow(I), viscircles(centers_desc_list, radii_desc_list,'Color','b'),rectangle('Position',[xmin_desc ymin_desc xmax_desc-xmin_desc ymax_desc-ymin_desc])
    
    %% Ascending Aorta
    %- Define rectangle
    y_width = size(I,1)/12.5; x_width = size(I,2)/12.5;
    xmin_asc = centers_asc_ref(1)-x_width; xmax_asc = centers_asc_ref(1)+x_width;
    ymin_asc = centers_asc_ref(2)-y_width; ymax_asc = centers_asc_ref(2)+x_width;
    %- Find potential descending aorta
    centers_asc_list = centers; radii_asc_list = radii;
    for k=size(centers_asc_list,1):-1:1
        if centers_asc_list(k,1)<xmin_asc || centers_asc_list(k,1)>xmax_asc || centers_asc_list(k,2)<ymin_asc || centers_asc_list(k,2)>ymax_asc
            centers_asc_list(k,:) = [];
            radii_asc_list(k) = [];
        end
    end
    %- Delete too large circles
    ind = find(radii_asc_list>15)
    radii_asc_list(ind) = [];
    centers_asc_list(ind,:) = [];
%     figure(),imshow(I), viscircles(centers_asc_list, radii_asc_list,'Color','b'),rectangle('Position',[xmin_asc ymin_asc xmax_asc-xmin_asc ymax_asc-ymin_asc])
    
    %% Plot
    if flag_plot ==1
        figure(),hold all,
        imshow(I),viscircles(centers_desc_list, radii_desc_list,'Color','b'),viscircles(centers_asc_list, radii_asc_list,'Color','r')
        rectangle('Position',[xmin_desc ymin_desc xmax_desc-xmin_desc ymax_desc-ymin_desc]),rectangle('Position',[xmin_asc ymin_asc xmax_asc-xmin_asc ymax_asc-ymin_asc])
    end
    
function [centers_desc radii_desc] = find_desc_from_ref(I,centers,centers_ref,radii,xmin_desc,xmax_desc,ymin_desc,ymax_desc,xmin_asc,flag_plot);
    %% Define middle location
    %- Y
    y_width = size(I,1)/7.5;
    ymid = size(I,1)/2;
    ylow = ymid*1.25 - y_width; yhigh = ymid*1.25 + y_width;
    %- X
    x_width = size(I,2)/7.5;
    xmid = size(I,2)/2;
    xlow = xmid - x_width; xhigh = xmid + x_width;
%     imshow(I),rectangle(end'Position',[xlow ylow xhigh-xlow yhigh-ylow])
    
    %% Select circles within the middle rectangle
    for k=size(centers,1):-1:1
        if centers(k,1)<xlow || centers(k,1)>xhigh || centers(k,2)<ylow || centers(k,2)>yhigh
            centers(k,:) = [];
            radii(k,:) = [];
        end
    end
    
    centers_desc_distance = distance_from_reference(centers,centers_ref);
    
    if flag_plot==1
        figure(),imshow(I), viscircles(centers, radii,'Color','b'),rectangle('Position',[xlow ylow xhigh-xlow yhigh-ylow])
    end
    
    %% Select descending aorta
    [c ind] = min(centers_desc_distance);
    centers_desc  = centers(ind,:);
    radii_desc = radii(ind);

function dataout = scaledata(datain,minval,maxval)
dataout = datain - min(datain(:));
dataout = (dataout/range(dataout(:)))*(maxval-minval);
dataout = dataout + minval;

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
    ind_separate = [0 ; ind_separate ; length(d_id_ok)];
    for k=1:2
        Time_ind{k} = [id_ok(ind_separate(k)+1):id_ok(ind_separate(k+1))]';
    end
else
    ind_separate = [0 ; ind_separate ; length(d_id_ok)];
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
elseif length(idx)==3
    variable_output_A = [{'Aasc'},{'Adesc'},{'Adia1'}];
    variable_output_Q = [{'Qasc'},{'Qdesc'},{'Qdia1'}];
    variable_output_t = [{'Tasc'},{'Tdesc'},{'Tdia1'}];
else
    variable_output_A = [{'Aasc'},{'Adesc'},{'Adia1'},{'Adia2'}];
    variable_output_Q = [{'Qasc'},{'Qdesc'},{'Qdia1'},{'Qdia2'}];
    variable_output_t = [{'Tasc'},{'Tdesc'},{'Tdia1'},{'Tdia2'}];
end

if length(idx)==2
    for k=1:length(idx)
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
else
    for k=1:length(idx)
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
    clear Adia1 Adia2 Qdia1 Qdia2 Tdia1 Tdia2 
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

function [] = generate_spreadsheet_1_part(filename,path,info,EF1,SV,HR,CO,length_asc_desc,transit_time_asc_desc,PWV_asc_desc,Dis_Asc,Dis_Desc,Tasc,Qasc,Qdesc,Aasc,Adesc,EDV,EF,PP,DBP,SBP,Pes,ENd_est,Ees,tNd,ENd_avg,C_total)

cd(path)

fill = {'','','','','',''};

%% PATIENT INFO
%- Check fields exist
if isfield(info,'PatientName')
    Field_Name = [info.PatientName.GivenName,' ',info.PatientName.FamilyName];
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
    Field_DOB = [info.PatientBirthDate(end-1:end),'-',info.PatientBirthDate(end-3:end-2),'-',info.PatientBirthDate(1:4)];
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

   %% FLOW ANALYSIS
line_Flow_Analysis = [{'%%% FLOW ANALYSIS','','','','',''};
       {'SV:',round(SV),'ml','','',''};
       {'HR:',round(HR),'bpm','','',''};
       {'CO:',round(100*CO)/100,'L/min','','',''};
       {'EDV:',round(EDV),'ml','','',''};
       {'EF:',round(EF),'%','','',''};
       {'EF1:',round(EF1),'%','','',''};
       ];
   
%% BLOOD PRESSURE
line_Pressure = [{'%%% BLOOD PRESSURE','','','','',''};
       {'DBP:',round(DBP),'mmHg','','',''};
       {'SBP:',round(SBP),'mmHg','','',''};
       {'PP:',round(PP),'mmHg','','',''};
       {'Pes:',round(Pes),'mmHg','','',''};
       ];   

%% ASCENDING-DESCENDING PWV
line_Asc_Desc_PWV = [{'%%% ASCENDING-DESCENDING PWV','','','','',''};
       {'Path length:',round(length_asc_desc)/1e3,'m','','',''};
       {'Transit time:',round(transit_time_asc_desc*1e4)/1e4,'s','','',''};
       {'PWV:',round(PWV_asc_desc*10)/10,'m/s','','',''};
       ];   
      
%% LV Elastance
line_Elastance = [{'%%% LEFT VENTRICULAR ELASTANCE','','','','',''};
       {'Elastance at onset of ejection (ENd):',round(ENd_est*100)/100,'mmHg/ml','','',''};
       {'Elastance at end systole (Ees):',round(Ees*100)/100,'mmHg/ml','','',''};
       ];   
   
%% Arterial Distensibility
line_Distensibility = [{'%%% ARTERIAL DISTENSIBILITY','','','','',''};
       {'Ascending aorta:',round(Dis_Asc*1e4)/1e4,'/mmHg','','',''};
       {'Descending aorta:',round(Dis_Desc*1e4)/1e4,'/mmHg','','',''};
       ];   
   
%% Flow waves
filler = strings(length(Qasc),1)
line_Flow_Waves = [["%%% FLOW WAVES","","","","",""];["Time (s)","Ascending aorta (ml/s)","Descending aorta (ml/s)","","",""];...
    [Tasc/1e3,Qasc,Qdesc,filler,filler,filler]]

%% Luminal variations waves
filler = strings(length(Aasc),1)
line_Luminal_Variations = [["%%% LUMINAL VARIATIONS","","","","",""];["Time (s)","Ascending aorta (mm2)","Descending aorta (mm2)","","",""];...
    [Tasc/1e3,Aasc,Adesc,filler,filler,filler]]

%% ASSEMBLE TABLE
T = [line_Patient_Info ; fill ; line_Flow_Analysis ; fill ; line_Pressure ; fill ; line_Asc_Desc_PWV ; fill ; line_Desc_Dia_PWV ; fill ; line_Asc_Dia_PWV ; fill ; line_Elastance ; fill ; line_Distensibility ; fill ; cellstr(line_Flow_Waves) ; fill ; cellstr(line_Luminal_Variations)];

%% WRITE FILE
writecell(T,[path filename]);

function [] = generate_spreadsheet_3_parts(filename,path,info,EF1,SV,HR,CO,length_asc_desc,length_desc_dia,length_asc_dia,transit_time_asc_desc,transit_time_desc_dia,transit_time_asc_dia,PWV_asc_desc,PWV_desc_dia,PWV_asc_dia,Dis_Asc,Dis_Desc,Dis_Dia,Tasc,Qasc,Qdesc,Qdia,Aasc,Adesc,Adia,EDV,EF,PP,DBP,SBP,Pes,ENd_est,Ees,tNd,ENd_avg,C_total)

cd(path)

fill = {'','','','','',''};

%% PATIENT INFO
%- Check fields exist
if isfield(info,'PatientName')
    Field_Name = [info.PatientName.GivenName,' ',info.PatientName.FamilyName];
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
    Field_DOB = [info.PatientBirthDate(end-1:end),'-',info.PatientBirthDate(end-3:end-2),'-',info.PatientBirthDate(1:4)];
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

   %% FLOW ANALYSIS
line_Flow_Analysis = [{'%%% FLOW ANALYSIS','','','','',''};
       {'SV:',round(SV),'ml','','',''};
       {'HR:',round(HR),'bpm','','',''};
       {'CO:',round(100*CO)/100,'L/min','','',''};
       {'EDV:',round(EDV),'ml','','',''};
       {'EF:',round(EF),'%','','',''};
       {'EF1:',round(EF1),'%','','',''};
       ];
   
%% BLOOD PRESSURE
line_Pressure = [{'%%% BLOOD PRESSURE','','','','',''};
       {'DBP:',round(DBP),'mmHg','','',''};
       {'SBP:',round(SBP),'mmHg','','',''};
       {'PP:',round(PP),'mmHg','','',''};
       {'Pes:',round(Pes),'mmHg','','',''};
       ];   

%% ASCENDING-DESCENDING PWV
line_Asc_Desc_PWV = [{'%%% ASCENDING-DESCENDING PWV','','','','',''};
       {'Path length:',round(length_asc_desc)/1e3,'m','','',''};
       {'Transit time:',round(transit_time_asc_desc*1e4)/1e4,'s','','',''};
       {'PWV:',round(PWV_asc_desc*10)/10,'m/s','','',''};
       ];   
   
%% DESCENDING-DIAPHRAGM PWV
line_Desc_Dia_PWV = [{'%%% DESCENDING-DIAPHRAGM PWV','','','','',''};
       {'Path length:',round(length_desc_dia)/1e3,'m','','',''};
       {'Transit time:',round(transit_time_desc_dia*1e4)/1e4,'s','','',''};
       {'PWV:',round(PWV_desc_dia*10)/10,'m/s','','',''};
       ];   
   
%% ASCENDING-DESCENDING PWV
line_Asc_Dia_PWV = [{'%%% ASCENDING-DIAPHRAGM PWV','','','','',''};
       {'Path length:',round(length_asc_dia)/1e3,'m','','',''};
       {'Transit time:',round(transit_time_asc_dia*1e4)/1e4,'s','','',''};
       {'PWV:',round(PWV_asc_dia*10)/10,'m/s','','',''};
       ];   
   
%% LV Elastance
line_Elastance = [{'%%% LEFT VENTRICULAR ELASTANCE','','','','',''};
       {'Elastance at onset of ejection (ENd):',round(ENd_est*100)/100,'mmHg/ml','','',''};
       {'Elastance at end systole (Ees):',round(Ees*100)/100,'mmHg/ml','','',''};
       ];   
   
%% Arterial Distensibility
line_Distensibility = [{'%%% ARTERIAL DISTENSIBILITY','','','','',''};
       {'Ascending aorta:',round(Dis_Asc*1e4)/1e4,'/mmHg','','',''};
       {'Descending aorta:',round(Dis_Desc*1e4)/1e4,'/mmHg','','',''};
       {'Diaphragm:',round(Dis_Dia*1e4)/1e4,'/mmHg','','',''};
       ];   
   
%% Flow waves
filler = strings(length(Qasc),1)
line_Flow_Waves = [["%%% FLOW WAVES","","","","",""];["Time (s)","Ascending aorta (ml/s)","Descending aorta (ml/s)","Diaphragm (ml/s)","",""];...
    [Tasc/1e3,Qasc,Qdesc,Qdia,filler,filler]]

%% Luminal variations waves
filler = strings(length(Aasc),1)
line_Luminal_Variations = [["%%% LUMINAL VARIATIONS","","","","",""];["Time (s)","Ascending aorta (mm2)","Descending aorta (mm2)","Diaphragm (mm2)","",""];...
    [Tasc/1e3,Aasc,Adesc,Adia,filler,filler]]

%% ASSEMBLE TABLE
T = [line_Patient_Info ; fill ; line_Flow_Analysis ; fill ; line_Pressure ; fill ; line_Asc_Desc_PWV ; fill ; line_Desc_Dia_PWV ; fill ; line_Asc_Dia_PWV ; fill ; line_Elastance ; fill ; line_Distensibility ; fill ; cellstr(line_Flow_Waves) ; fill ; cellstr(line_Luminal_Variations)];

%% WRITE FILE
writecell(T,[path filename]);


function ObjectPolarity_mode = decide_ObjectPolarity_mode(im);
    m = mean(mean(im));
    if m <0.35
        ObjectPolarity_mode = 'bright';
    else
        ObjectPolarity_mode = 'dark';
    end
    
function [pt,dudt,fofthandle] = interparc(t,px,py,varargin)
% interparc: interpolate points along a curve in 2 or more dimensions
% usage: pt = interparc(t,px,py)    % a 2-d curve
% usage: pt = interparc(t,px,py,pz) % a 3-d curve
% usage: pt = interparc(t,px,py,pz,pw,...) % a 4-d or higher dimensional curve
% usage: pt = interparc(t,px,py,method) % a 2-d curve, method is specified
% usage: [pt,dudt,fofthandle] = interparc(t,px,py,...) % also returns derivatives, and a function handle
%
% Interpolates new points at any fractional point along
% the curve defined by a list of points in 2 or more
% dimensions. The curve may be defined by any sequence
% of non-replicated points.
%
% arguments: (input)
%  t   - vector of numbers, 0 <= t <= 1, that define
%        the fractional distance along the curve to
%        interpolate the curve at. t = 0 will generate
%        the very first point in the point list, and
%        t = 1 yields the last point in that list.
%        Similarly, t = 0.5 will yield the mid-point
%        on the curve in terms of arc length as the
%        curve is interpolated by a parametric spline.
%
%        If t is a scalar integer, at least 2, then
%        it specifies the number of equally spaced
%        points in arclength to be generated along
%        the curve.
%
%  px, py, pz, ... - vectors of length n, defining
%        points along the curve. n must be at least 2.
%        Exact Replicate points should not be present
%        in the curve, although there is no constraint
%        that the curve has replicate independent
%        variables.
%
%  method - (OPTIONAL) string flag - denotes the method
%        used to compute the points along the curve.
%
%        method may be any of 'linear', 'spline', or 'pchip',
%        or any simple contraction thereof, such as 'lin',
%        'sp', or even 'p'.
%        
%        method == 'linear' --> Uses a linear chordal
%               approximation to interpolate the curve.
%               This method is the most efficient.
%
%        method == 'pchip' --> Uses a parametric pchip
%               approximation for the interpolation
%               in arc length.
%
%        method == 'spline' --> Uses a parametric spline
%               approximation for the interpolation in
%               arc length. Generally for a smooth curve,
%               this method may be most accurate.
%
%        method = 'csape' --> if available, this tool will
%               allow a periodic spline fit for closed curves.
%               ONLY use this method if your points should
%               represent a closed curve.
%               
%               If the last point is NOT the same as the
%               first point on the curve, then the curve
%               will be forced to be periodic by this option.
%               That is, the first point will be replicated
%               onto the end.
%
%               If csape is not present in your matlab release,
%               then an error will result.
%
%        DEFAULT: 'spline'
%
%
% arguments: (output)
%  pt - Interpolated points at the specified fractional
%        distance (in arc length) along the curve.
%
%  dudt - when a second return argument is required,
%       interparc will return the parametric derivatives
%       (dx/dt, dy/dt, dz/dt, ...) as an array.
%
%  fofthandle - a function handle, taking numbers in the interval [0,1]
%       and evaluating the function at those points.
%
%       Extrapolation will not be permitted by this call.
%       Any values of t that lie outside of the interval [0,1]
%       will be clipped to the endpoints of the curve.
%
% Example:
% % Interpolate a set of unequally spaced points around
% % the perimeter of a unit circle, generating equally
% % spaced points around the perimeter.
% theta = sort(rand(15,1))*2*pi;
% theta(end+1) = theta(1);
% px = cos(theta);
% py = sin(theta);
%
% % interpolate using parametric splines
% pt = interparc(100,px,py,'spline');
%
% % Plot the result
% plot(px,py,'r*',pt(:,1),pt(:,2),'b-o')
% axis([-1.1 1.1 -1.1 1.1])
% axis equal
% grid on
% xlabel X
% ylabel Y
% title 'Points in blue are uniform in arclength around the circle'
%
%
% Example:
% % For the previous set of points, generate exactly 6
% % points around the parametric splines, verifying
% % the uniformity of the arc length interpolant.
% pt = interparc(6,px,py,'spline');
%
% % Convert back to polar form. See that the radius
% % is indeed 1, quite accurately.
% [TH,R] = cart2pol(pt(:,1),pt(:,2))
% % TH =
% %       0.86005
% %        2.1141
% %       -2.9117
% %        -1.654
% %      -0.39649
% %       0.86005
% % R =
% %             1
% %        0.9997
% %        0.9998
% %       0.99999
% %        1.0001
% %             1
%
% % Unwrap the polar angles, and difference them.
% diff(unwrap(TH))
% % ans =
% %        1.2541
% %        1.2573
% %        1.2577
% %        1.2575
% %        1.2565
%
% % Six points around the circle should be separated by
% % 2*pi/5 radians, if they were perfectly uniform. The
% % slight differences are due to the imperfect accuracy
% % of the parametric splines.
% 2*pi/5
% % ans =
% %        1.2566
%
%
% See also: arclength, spline, pchip, interp1
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 3/15/2010

% unpack the arguments and check for errors
if nargin < 3
  error('ARCLENGTH:insufficientarguments', ...
    'at least t, px, and py must be supplied')
end

t = t(:);
if (numel(t) == 1) && (t > 1) && (rem(t,1) == 0)
  % t specifies the number of points to be generated
  % equally spaced in arclength
  t = linspace(0,1,t)';
elseif any(t < 0) || any(t > 1)
  error('ARCLENGTH:impropert', ...
    'All elements of t must be 0 <= t <= 1')
end

% how many points will be interpolated?
nt = numel(t);

% the number of points on the curve itself
px = px(:);
py = py(:);
n = numel(px);

% are px and py both vectors of the same length?
if ~isvector(px) || ~isvector(py) || (length(py) ~= n)
  error('ARCLENGTH:improperpxorpy', ...
    'px and py must be vectors of the same length')
elseif n < 2
  error('ARCLENGTH:improperpxorpy', ...
    'px and py must be vectors of length at least 2')
end

% compose px and py into a single array. this way,
% if more dimensions are provided, the extension
% is trivial.
pxy = [px,py];
ndim = 2;

% the default method is 'linear'
method = 'spline';

% are there any other arguments?
if nargin > 3
  % there are. check the last argument. Is it a string?
  if ischar(varargin{end})
    method = varargin{end};
    varargin(end) = [];
    
    % method may be any of {'linear', 'pchip', 'spline', 'csape'.}
    % any any simple contraction thereof.
    valid = {'linear', 'pchip', 'spline', 'csape'};
    [method,errstr] = validstring(method,valid);
    if ~isempty(errstr)
      error('INTERPARC:incorrectmethod',errstr)
    end
  end
  
  % anything that remains in varargin must add
  % an additional dimension on the curve/polygon
  for i = 1:numel(varargin)
    pz = varargin{i};
    pz = pz(:);
    if numel(pz) ~= n
      error('ARCLENGTH:improperpxorpy', ...
        'pz must be of the same size as px and py')
    end
    pxy = [pxy,pz]; %#ok
  end
  
  % the final number of dimensions provided
  ndim = size(pxy,2);
end

% if csape, then make sure the first point is replicated at the end.
% also test to see if csape is available
if method(1) == 'c'
  if exist('csape','file') == 0
    error('CSAPE was requested, but you lack the necessary toolbox.')
  end
  
  p1 = pxy(1,:);
  pend = pxy(end,:);
  
  % get a tolerance on whether the first point is replicated.
  if norm(p1 - pend) > 10*eps(norm(max(abs(pxy),[],1)))
    % the two end points were not identical, so wrap the curve
    pxy(end+1,:) = p1;
    nt = nt + 1;
  end
end

% preallocate the result, pt
pt = NaN(nt,ndim);

% Compute the chordal (linear) arclength
% of each segment. This will be needed for
% any of the methods.
chordlen = sqrt(sum(diff(pxy,[],1).^2,2));

% Normalize the arclengths to a unit total
chordlen = chordlen/sum(chordlen);

% cumulative arclength
cumarc = [0;cumsum(chordlen)];

% The linear interpolant is trivial. do it as a special case
if method(1) == 'l'
  % The linear method.
  
  % which interval did each point fall in, in
  % terms of t?
  [junk,tbins] = histc(t,cumarc); %#ok
  
  % catch any problems at the ends
  tbins((tbins <= 0) | (t <= 0)) = 1;
  tbins((tbins >= n) | (t >= 1)) = n - 1;
  
  % interpolate
  s = (t - cumarc(tbins))./chordlen(tbins);
  % be nice, and allow the code to work on older releases
  % that don't have bsxfun
  pt = pxy(tbins,:) + (pxy(tbins+1,:) - pxy(tbins,:)).*repmat(s,1,ndim);
  
  % do we need to compute derivatives here?
  if nargout > 1
    dudt = (pxy(tbins+1,:) - pxy(tbins,:))./repmat(chordlen(tbins),1,ndim);
  end
  
  % do we need to create the spline as a piecewise linear function?
  if nargout > 2
    spl = cell(1,ndim);
    for i = 1:ndim
      coefs = [diff(pxy(:,i))./diff(cumarc),pxy(1:(end-1),i)];
      spl{i} = mkpp(cumarc.',coefs);
    end
    
    %create a function handle for evaluation, passing in the splines
    fofthandle = @(t) foft(t,spl);
  end
  
  % we are done at this point
  return
end

% If we drop down to here, we have either a spline
% or csape or pchip interpolant to work with.

% compute parametric splines
spl = cell(1,ndim);
spld = spl;
diffarray = [3 0 0;0 2 0;0 0 1;0 0 0];
for i = 1:ndim
  switch method
    case 'pchip'
      spl{i} = pchip(cumarc,pxy(:,i));
    case 'spline'
      spl{i} = spline(cumarc,pxy(:,i));
      nc = numel(spl{i}.coefs);
      if nc < 4
        % just pretend it has cubic segments
        spl{i}.coefs = [zeros(1,4-nc),spl{i}.coefs];
        spl{i}.order = 4;
      end
    case 'csape'
      % csape was specified, so the curve is presumed closed,
      % therefore periodic
      spl{i} = csape(cumarc,pxy(:,i),'periodic');
      nc = numel(spl{i}.coefs);
      if nc < 4
        % just pretend it has cubic segments
        spl{i}.coefs = [zeros(1,4-nc),spl{i}.coefs];
        spl{i}.order = 4;
      end
  end
  
  % and now differentiate them
  xp = spl{i};
  xp.coefs = xp.coefs*diffarray;
  xp.order = 3;
  spld{i} = xp;
end

% catch the case where there were exactly three points
% in the curve, and spline was used to generate the
% interpolant. In this case, spline creates a curve with
% only one piece, not two.
if (numel(cumarc) == 3) && (method(1) == 's')
  cumarc = spl{1}.breaks;
  n = numel(cumarc);
  chordlen = sum(chordlen);
end

% Generate the total arclength along the curve
% by integrating each segment and summing the
% results. The integration scheme does its job
% using an ode solver.

% polyarray here contains the derivative polynomials
% for each spline in a given segment
polyarray = zeros(ndim,3);
seglen = zeros(n-1,1);

% options for ode45
opts = odeset('reltol',1.e-9);
for i = 1:spl{1}.pieces
  % extract polynomials for the derivatives
  for j = 1:ndim
    polyarray(j,:) = spld{j}.coefs(i,:);
  end
  
  % integrate the arclength for the i'th segment
  % using ode45 for the integral. I could have
  % done this part with quad too, but then it
  % would not have been perfectly (numerically)
  % consistent with the next operation in this tool.
  [tout,yout] = ode45(@(t,y) segkernel(t,y,ndim,polyarray),[0,chordlen(i)],0,opts); %#ok
  seglen(i) = yout(end);
end

% and normalize the segments to have unit total length
totalsplinelength = sum(seglen);
cumseglen = [0;cumsum(seglen)];

% which interval did each point fall into, in
% terms of t, but relative to the cumulative
% arc lengths along the parametric spline?
[junk,tbins] = histc(t*totalsplinelength,cumseglen); %#ok

% catch any problems at the ends
tbins((tbins <= 0) | (t <= 0)) = 1;
tbins((tbins >= n) | (t >= 1)) = n - 1;

% Do the fractional integration within each segment
% for the interpolated points. t is the parameter
% used to define the splines. It is defined in terms
% of a linear chordal arclength. This works nicely when
% a linear piecewise interpolant was used. However,
% what is asked for is an arclength interpolation
% in terms of arclength of the spline itself. Call s
% the arclength traveled along the spline.
s = totalsplinelength*t;

% the ode45 options will now include an events property
% so we can catch zero crossings.
opts = odeset('reltol',1.e-9,'events',@ode_events);

ti = t;
for i = 1:nt
  % si is the piece of arc length that we will look
  % for in this spline segment.
  si = s(i) - cumseglen(tbins(i));
  
  % extract polynomials for the derivatives
  % in the interval the point lies in
  for j = 1:ndim
    polyarray(j,:) = spld{j}.coefs(tbins(i),:);
  end
  
  % we need to integrate in t, until the integral
  % crosses the specified value of si. Because we
  % have defined totalsplinelength, the lengths will
  % be normalized at this point to a unit length.
  %
  % Start the ode solver at -si, so we will just
  % look for an event where y crosses zero.
  [tout,yout,te,ye] = ode45(@(t,y) segkernel(t,y,ndim,polyarray),[0,chordlen(tbins(i))],-si,opts); %#ok
  
  % we only need that point where a zero crossing occurred
  % if no crossing was found, then we can look at each end.
  if ~isempty(te)
    ti(i) = te(1) + cumarc(tbins(i));
  else
    % a crossing must have happened at the very
    % beginning or the end, and the ode solver
    % missed it, not trapping that event.
    if abs(yout(1)) < abs(yout(end))
      % the event must have been at the start.
      ti(i) = tout(1) + cumarc(tbins(i));
    else
      % the event must have been at the end.
      ti(i) = tout(end) + cumarc(tbins(i));
    end
  end
end

% Interpolate the parametric splines at ti to get
% our interpolated value.
for L = 1:ndim
  pt(:,L) = ppval(spl{L},ti);
end

% do we need to compute first derivatives here at each point?
if nargout > 1
  dudt = zeros(nt,ndim);
  for L = 1:ndim
    dudt(:,L) = ppval(spld{L},ti);
  end
end

% create a function handle for evaluation, passing in the splines
if nargout > 2
  fofthandle = @(t) foft(t,spl);
end

function val = segkernel(t,y,ndim,polyarray) %#ok
    % sqrt((dx/dt)^2 + (dy/dt)^2 + ...)
    val = zeros(size(t));
    for k = 1:ndim
      val = val + polyval(polyarray(k,:),t).^2;
    end
    val = sqrt(val);

function [value,isterminal,direction] = ode_events(t,y) %#ok
  % ode event trap, looking for zero crossings of y.
  value = y;
  isterminal = ones(size(y));
  direction = ones(size(y));

function f_t = foft(t,spl)
% tool allowing the user to evaluate the interpolant at any given point for any values t in [0,1]
pdim = numel(spl);
f_t = zeros(numel(t),pdim);

% convert t to a column vector, clipping it to [0,1] as we do.
t = max(0,min(1,t(:)));

% just loop over the splines in the cell array of splines
for i = 1:pdim
  f_t(:,i) = ppval(spl{i},t);
end

function [str,errorclass] = validstring(arg,valid)
% validstring: compares a string against a set of valid options
% usage: [str,errorclass] = validstring(arg,valid)
%
% If a direct hit, or any unambiguous shortening is found, that
% string is returned. Capitalization is ignored.
%
% arguments: (input)
%  arg - character string, to be tested against a list
%        of valid choices. Capitalization is ignored.
%
%  valid - cellstring array of alternative choices
%
% Arguments: (output)
%  str - string - resulting choice resolved from the
%        list of valid arguments. If no unambiguous
%        choice can be resolved, then str will be empty.
%
%  errorclass - string - A string argument that explains
%        the error. It will be one of the following
%        possibilities:
%
%        ''  --> No error. An unambiguous match for arg
%                was found among the choices.
%
%        'No match found' --> No match was found among 
%                the choices provided in valid.
%
%        'Ambiguous argument' --> At least two ambiguous
%                matches were found among those provided
%                in valid.
%        
%
% Example:
%  valid = {'off' 'on' 'The sky is falling'}
%  
%
% See also: parse_pv_pairs, strmatch, strcmpi
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 3/25/2010

ind = find(strncmpi(lower(arg),valid,numel(arg)));
if isempty(ind)
  % No hit found
  errorclass = 'No match found';
  str = '';
elseif (length(ind) > 1)
  % Ambiguous arg, hitting more than one of the valid options
  errorclass = 'Ambiguous argument';
  str = '';
  return
else
  errorclass = '';
  str = valid{ind};
end

