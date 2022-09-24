function varargout = GUI_main(varargin)
% GUI_MAIN MATLAB code for GUI_main.fig
%      GUI_MAIN, by itself, creates a new GUI_MAIN or raises the existing
%      singleton*.
%
%      H = GUI_MAIN returns the handle to a new GUI_MAIN or the handle to
%      the existing singleton*.
%
%      GUI_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_MAIN.M with the given input arguments.
%
%      GUI_MAIN('Property','Value',...) creates a new GUI_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_main

% Last Modified by GUIDE v2.5 22-Sep-2022 18:28:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_main_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_main_OutputFcn, ...
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


% --- Executes just before GUI_main is made visible.
function GUI_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_main (see VARARGIN)

% Choose default command line output for GUI_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_axial.
function button_axial_Callback(hObject, eventdata, handles)
% hObject    handle to button_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_pick_axial_mode

% --- Executes on button press in button_parasagittal.
function button_parasagittal_Callback(hObject, eventdata, handles)
% hObject    handle to button_parasagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_PWV_parasagittal_final

% --- Executes on button press in button_fuse.
% function button_fuse_Callback(hObject, eventdata, handles)
% % hObject    handle to button_fuse (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% GUI_PWV_parasagittal_final_fuse

% --- Executes on button press in button_3D.
function button_3D_Callback(hObject, eventdata, handles)
% hObject    handle to button_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_pick_coronal_mode


% --- Executes on button press in pushbutton_LA_LV_coupling.
function pushbutton_LA_LV_coupling_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_LA_LV_coupling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GUI_coupling_LA_LV
