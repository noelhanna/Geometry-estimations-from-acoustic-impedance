function varargout = GUI_Fit_v0_2(varargin)
% GUI_FIT_V0_2 MATLAB code for GUI_Fit_v0_2.fig
%      GUI_FIT_V0_2, by itself, creates a new GUI_FIT_V0_2 or raises the existing
%      singleton*.
%
%      H = GUI_FIT_V0_2 returns the handle to a new GUI_FIT_V0_2 or the handle to
%      the existing singleton*.
%
%      GUI_FIT_V0_2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FIT_V0_2.M with the given input arguments.
%
%      GUI_FIT_V0_2('Property','Value',...) creates a new GUI_FIT_V0_2 or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_Fit_v0_2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_Fit_v0_2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Fit_v0_2

% Last Modified by GUIDE v2.5 21-Nov-2018 13:38:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Fit_v0_2_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Fit_v0_2_OutputFcn, ...
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
% -----------------------------------------------


% --- Executes just before GUI_Fit_v0_2 is made visible.
function GUI_Fit_v0_2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Fit_v0_2 (see VARARGIN)

% Choose default command line output for GUI_Fit_v0_2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes GUI_Fit_v0_2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Fit_v0_2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end


% Update handles structure
guidata(handles.figure1, handles);

%% ---------------------------------------------------------

%% Vocal Tract

% --- Executes on selection change in VTsegments. Select number of VT
% cylinders
function VTsegments_Callback(hObject, eventdata, handles)
% hObject    handle to VTsegments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns VTsegments contents as cell array
%        contents{get(hObject,'Value')} returns selected item from VTsegments
contents = cellstr(get(hObject,'String'));
VT.segments = contents{get(hObject,'Value')};
handles.VTinitialLength.Value = 170/eval(VT.segments); % initial length assumes 170mm length
set(handles.VTinitialLength, 'String', num2str(handles.VTinitialLength.Value));


% --- Executes during object creation, after setting all properties.
function VTsegments_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VTsegments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function VTinitialRadius_Callback(hObject, eventdata, handles)
% hObject    handle to VTinitialRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VTinitialRadius as text
%        str2double(get(hObject,'String')) returns contents of VTinitialRadius as a double
contents = cellstr(get(hObject,'String'));
VT.radius = set(hObject,'Value', str2double(contents));


% --- Executes during object creation, after setting all properties.
function VTinitialRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VTinitialRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function VTminRadius_Callback(hObject, eventdata, handles)
% hObject    handle to VTminRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VTminRadius as text
%        str2double(get(hObject,'String')) returns contents of VTminRadius as a double


% --- Executes during object creation, after setting all properties.
function VTminRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VTminRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function VTmaxRadius_Callback(hObject, eventdata, handles)
% hObject    handle to VTmaxRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VTmaxRadius as text
%        str2double(get(hObject,'String')) returns contents of VTmaxRadius as a double


% --- Executes during object creation, after setting all properties.
function VTmaxRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VTmaxRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function VTinitialLength_Callback(hObject, eventdata, handles)
% hObject    handle to VTinitialLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VTinitialLength as text
%        str2double(get(hObject,'String')) returns contents of VTinitialLength as a double


% --- Executes during object creation, after setting all properties.
function VTinitialLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VTinitialLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% set(handles.VTinitialLength, 'String', num2str(handles.VTinitialLength.Value));


function VTminLength_Callback(hObject, eventdata, handles)
% hObject    handle to VTminLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VTminLength as text
%        str2double(get(hObject,'String')) returns contents of VTminLength as a double


% --- Executes during object creation, after setting all properties.
function VTminLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VTminLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function VTmaxLength_Callback(hObject, eventdata, handles)
% hObject    handle to VTmaxLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VTmaxLength as text
%        str2double(get(hObject,'String')) returns contents of VTmaxLength as a double


% --- Executes during object creation, after setting all properties.
function VTmaxLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VTmaxLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%% Glottis
% --- Executes on selection change in GlottisState.
function GlottisState_Callback(hObject, eventdata, handles)
% hObject    handle to GlottisState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns GlottisState contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GlottisState

if hObject.Value == 1 % glottis closed
    set(handles.GlottisInitialRadius, 'Value', 0.000001);
    set(handles.GlottisInitialRadius, 'String', '0.000001');
    
    set(handles.GlottisMinRadius, 'Value', -10);
    set(handles.GlottisMinRadius, 'String', '-10');
    
    set(handles.GlottisMaxRadius, 'Value', 0.000009);
    set(handles.GlottisMaxRadius, 'String', '0.000009');
        
else % glottis open
    set(handles.GlottisInitialRadius, 'Value', 1);
    set(handles.GlottisInitialRadius, 'String', '1');
    
    set(handles.GlottisMinRadius, 'Value', 0.01);
    set(handles.GlottisMinRadius, 'String', '0.01');
    
    set(handles.GlottisMaxRadius, 'Value', 20);
    set(handles.GlottisMaxRadius, 'String', '20');
    
    set(handles.VTsegments, 'Value', 7);

end




% --- Executes during object creation, after setting all properties.
function GlottisState_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GlottisState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function GlottisMaxRadius_Callback(hObject, eventdata, handles)
% hObject    handle to GlottisMaxRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GlottisMaxRadius as text
%        str2double(get(hObject,'String')) returns contents of GlottisMaxRadius as a double


% --- Executes during object creation, after setting all properties.
function GlottisMaxRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GlottisMaxRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function GlottisMaxLength_Callback(hObject, eventdata, handles)
% hObject    handle to GlottisMaxLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GlottisMaxLength as text
%        str2double(get(hObject,'String')) returns contents of GlottisMaxLength as a double


% --- Executes during object creation, after setting all properties.
function GlottisMaxLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GlottisMaxLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function GlottisMinLength_Callback(hObject, eventdata, handles)
% hObject    handle to GlottisMinLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GlottisMinLength as text
%        str2double(get(hObject,'String')) returns contents of GlottisMinLength as a double


% --- Executes during object creation, after setting all properties.
function GlottisMinLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GlottisMinLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function GlottisInitialLength_Callback(hObject, eventdata, handles)
% hObject    handle to GlottisInitialLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GlottisInitialLength as text
%        str2double(get(hObject,'String')) returns contents of GlottisInitialLength as a double


% --- Executes during object creation, after setting all properties.
function GlottisInitialLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GlottisInitialLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function GlottisMinRadius_Callback(hObject, eventdata, handles)
% hObject    handle to GlottisMinRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GlottisMinRadius as text
%        str2double(get(hObject,'String')) returns contents of GlottisMinRadius as a double


% --- Executes during object creation, after setting all properties.
function GlottisMinRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GlottisMinRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function GlottisInitialRadius_Callback(hObject, eventdata, handles)
% hObject    handle to GlottisInitialRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GlottisInitialRadius as text
%        str2double(get(hObject,'String')) returns contents of GlottisInitialRadius as a double


% --- Executes during object creation, after setting all properties.
function GlottisInitialRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GlottisInitialRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%% Subglottal
function SGinitialRadius_Callback(hObject, eventdata, handles)
% hObject    handle to SGinitialRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SGinitialRadius as text
%        str2double(get(hObject,'String')) returns contents of SGinitialRadius as a double


% --- Executes during object creation, after setting all properties.
function SGinitialRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SGinitialRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SGminRadius_Callback(hObject, eventdata, handles)
% hObject    handle to SGminRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SGminRadius as text
%        str2double(get(hObject,'String')) returns contents of SGminRadius as a double


% --- Executes during object creation, after setting all properties.
function SGminRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SGminRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SGinitialLength_Callback(hObject, eventdata, handles)
% hObject    handle to SGinitialLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SGinitialLength as text
%        str2double(get(hObject,'String')) returns contents of SGinitialLength as a double


% --- Executes during object creation, after setting all properties.
function SGinitialLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SGinitialLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SGminLength_Callback(hObject, eventdata, handles)
% hObject    handle to SGminLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SGminLength as text
%        str2double(get(hObject,'String')) returns contents of SGminLength as a double


% --- Executes during object creation, after setting all properties.
function SGminLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SGminLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SGmaxLength_Callback(hObject, eventdata, handles)
% hObject    handle to SGmaxLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SGmaxLength as text
%        str2double(get(hObject,'String')) returns contents of SGmaxLength as a double


% --- Executes during object creation, after setting all properties.
function SGmaxLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SGmaxLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SGmaxRadius_Callback(hObject, eventdata, handles)
% hObject    handle to SGmaxRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SGmaxRadius as text
%        str2double(get(hObject,'String')) returns contents of SGmaxRadius as a double


% --- Executes during object creation, after setting all properties.
function SGmaxRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SGmaxRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SGsegments.
function SGsegments_Callback(hObject, eventdata, handles)
% hObject    handle to SGsegments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SGsegments contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SGsegments
contents = cellstr(get(hObject,'String'));
SGsegments = contents{get(hObject,'Value')};
handles.SGinitialLength.Value = 200/eval(SGsegments); % initial length assumes 20cm length
set(handles.SGinitialLength, 'String', num2str(handles.SGinitialLength.Value));


% --- Executes during object creation, after setting all properties.
function SGsegments_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SGsegments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Non-Rigid 
% --- Executes on selection change in Rigidity.
function Rigidity_Callback(hObject, eventdata, handles)
% hObject    handle to Rigidity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Rigidity contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Rigidity


% --- Executes during object creation, after setting all properties.
function Rigidity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rigidity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function LmaxInertance_Callback(hObject, eventdata, handles)
% hObject    handle to LmaxInertance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LmaxInertance as text
%        str2double(get(hObject,'String')) returns contents of LmaxInertance as a double


% --- Executes during object creation, after setting all properties.
function LmaxInertance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LmaxInertance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function AmaxAttenuation_Callback(hObject, eventdata, handles)
% hObject    handle to AmaxAttenuation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AmaxAttenuation as text
%        str2double(get(hObject,'String')) returns contents of AmaxAttenuation as a double


% --- Executes during object creation, after setting all properties.
function AmaxAttenuation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AmaxAttenuation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function AminAttenuation_Callback(hObject, eventdata, handles)
% hObject    handle to AminAttenuation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AminAttenuation as text
%        str2double(get(hObject,'String')) returns contents of AminAttenuation as a double


% --- Executes during object creation, after setting all properties.
function AminAttenuation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AminAttenuation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function AinitialAttenuation_Callback(hObject, eventdata, handles)
% hObject    handle to AinitialAttenuation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AinitialAttenuation as text
%        str2double(get(hObject,'String')) returns contents of AinitialAttenuation as a double


% --- Executes during object creation, after setting all properties.
function AinitialAttenuation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AinitialAttenuation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LminIntertance_Callback(hObject, eventdata, handles)
% hObject    handle to LminIntertance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LminIntertance as text
%        str2double(get(hObject,'String')) returns contents of LminIntertance as a double


% --- Executes during object creation, after setting all properties.
function LminIntertance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LminIntertance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function LinitialInertance_Callback(hObject, eventdata, handles)
% hObject    handle to LinitialInertance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LinitialInertance as text
%        str2double(get(hObject,'String')) returns contents of LinitialInertance as a double


% --- Executes during object creation, after setting all properties.
function LinitialInertance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LinitialInertance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% update all values to match their edited strings
handles.VTinitialRadius.Value = str2double(handles.VTinitialRadius.String);
handles.VTmaxRadius.Value = str2double(handles.VTmaxRadius.String);
handles.VTminRadius.Value = str2double(handles.VTminRadius.String);
handles.VTinitialLength.Value = str2double(handles.VTinitialLength.String);
handles.VTmaxLength.Value = str2double(handles.VTmaxLength.String);
handles.VTminLength.Value = str2double(handles.VTminLength.String);
handles.SGinitialRadius.Value = str2double(handles.SGinitialRadius.String);
handles.SGmaxRadius.Value = str2double(handles.SGmaxRadius.String);
handles.SGminRadius.Value = str2double(handles.SGminRadius.String);
handles.SGinitialLength.Value = str2double(handles.SGinitialLength.String);
handles.SGmaxLength.Value = str2double(handles.SGmaxLength.String);
handles.SGminLength.Value = str2double(handles.SGminLength.String);
handles.GlottisInitialRadius.Value = str2double(handles.GlottisInitialRadius.String);
handles.GlottisMaxRadius.Value = str2double(handles.GlottisMaxRadius.String);
handles.GlottisMinRadius.Value = str2double(handles.GlottisMinRadius.String);
handles.GlottisInitialLength.Value = str2double(handles.GlottisInitialLength.String);
handles.GlottisMaxLength.Value = str2double(handles.GlottisMaxLength.String);
handles.GlottisMinLength.Value = str2double(handles.GlottisMinLength.String);
handles.AinitialAttenuation.Value = str2double(handles.AinitialAttenuation.String);
handles.AmaxAttenuation.Value = str2double(handles.AmaxAttenuation.String);
handles.AminAttenuation.Value = str2double(handles.AminAttenuation.String);




VT.segments = handles.VTsegments.Value;
VT.radius = handles.VTinitialRadius.Value.*ones(VT.segments,1)./1000; % m
VT.radiusLimits = [handles.VTminRadius.Value, handles.VTmaxRadius.Value]./1000;
VT.length = handles.VTinitialLength.Value.*ones(VT.segments,1)./1000;
VT.lengthLimits = [handles.VTminLength.Value, handles.VTmaxLength.Value]./1000;
VT.alpha_multiplier = handles.AinitialAttenuation.Value.*ones(VT.segments,1);

glottis.state = handles.GlottisState.Value;
glottis.radius = handles.GlottisInitialRadius.Value./1000;
glottis.length = handles.GlottisInitialLength.Value./1000;
glottis.radiusLimits = [handles.GlottisMinRadius.Value, handles.GlottisMaxRadius.Value]./1000;
glottis.lengthLimits = [handles.GlottisMinLength.Value, handles.GlottisMaxLength.Value]./1000;
glottis.alpha_multiplier = handles.AinitialAttenuation.Value;

SG.segments = handles.SGsegments.Value;
SG.radius = handles.SGinitialRadius.Value.*ones(SG.segments,1)./1000;
SG.radiusLimits = [handles.SGminRadius.Value, handles.SGmaxRadius.Value]./1000;
SG.length = handles.SGinitialLength.Value.*ones(SG.segments,1)./1000;
SG.lengthLimits = [handles.SGminLength.Value, handles.SGmaxLength.Value]./1000;
SG.alpha_multiplier = handles.AinitialAttenuation.Value.*ones(SG.segments,1);

Rigidity.rigidity = handles.Rigidity.Value;
Rigidity.inertanceInitial = handles.LinitialInertance.Value;
Rigidity.inertanceLimits = [handles.LminIntertance.Value, handles.LmaxInertance.Value];
Rigidity.alpha_multiplier = handles.AinitialAttenuation.Value;
Rigidity.alphaLimits = [handles.AminAttenuation.Value, handles.AmaxAttenuation.Value];



save('FitParameters', 'VT', 'glottis', 'SG', 'Rigidity');

fit_general_cylinder_model_v0_3(handles,...
    VT, glottis, SG, Rigidity);


% --- Executes on button press in LoadImpedanceFile.
function LoadImpedanceFile_Callback(hObject, eventdata, handles)
% hObject    handle to LoadImpedanceFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname, filterindex] = uigetfile('*.mat',...
    'Choose an impedance .mat file to fit:');

if filename == 0
    uiwait(msgbox('No file selected. Default used.', 'errorbox', 'error'));
    load('Closed_In_Ex.mat');
    filename = 'Closed_In_Ex.mat';
else
    load(filename);
end

axes(handles.axes2); hold off;
axes(handles.axes3); hold off;


set(handles.LoadImpedanceFile, 'String', filename);
set(handles.IterationSelection, 'Value', 1);

% load(filename)
% % load('Closed_In_Ex.mat')
noIterations = length(Iteration.Z(1,:));
set(handles.IterationSelection, 'String',...
    cellstr(string(1:1:noIterations)'));

Z_measured = Iteration.Z(:,handles.IterationSelection.Value);
x = Parameters.frequencyVector;

% get(handles.axes2);
plot(handles.axes2, x, abs(Z_measured), '.b'); 
set(handles.axes2, 'YScale', 'log');
xlabel(handles.axes2, 'Frequency (Hz)');
ylabel(handles.axes2, 'Z at the lips (10^{})');
axis tight

% get(handles.axes3);
plot(handles.axes3, x, angle(Z_measured), '.b'); 
set(handles.axes3 , 'YScale', 'linear');
xlabel(handles.axes3, 'Frequency (Hz)');
ylabel(handles.axes3, 'Angle (rad)');
ylim(handles.axes3, [-2 2]);
axis tight

linkaxes([handles.axes2, handles.axes3], 'x')


% --- Executes on selection change in IterationSelection.
function IterationSelection_Callback(hObject, eventdata, handles)
% hObject    handle to IterationSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns IterationSelection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from IterationSelection

load(handles.LoadImpedanceFile.String);
Z_measured = Iteration.Z(:,hObject.Value);
x = Parameters.frequencyVector;

axes(handles.axes2); hold off;
axes(handles.axes3); hold off;

plot(handles.axes2, x, abs(log10(Z_measured)), '.b');

plot(handles.axes3, x, angle(Z_measured), '.b'); 
ylim(handles.axes3, [-2 2]);


% --- Executes during object creation, after setting all properties.
function IterationSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IterationSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2
