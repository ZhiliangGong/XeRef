function varargout = molviewer_gui(mmsource, varargin)
%MOLVIEWER visualize molecules.
%
%   H = MOLVIEWER is an interactive tool for viewing a molecule 3D
%   structure.
%
%   H = MOLVIEWER(FILENAME) loads a molecule model file into the viewer,
%   and returns a handle H of the figure window that contains the viewer.
%   Supported file formats are PDB, CIF/mmCIF, MDL's mol, SDF, XYZ, SMOL
%   and JVXL, etc.
%
%   Note: Right click on the viewer to display commands for manipulating
%   the molecule display.
%
%   MOLVIEWER(PDBID) search for the PDBID in the Protein Data Bank and
%   displays the returned PDB model.
%
%   MOLVIEWER(URL) search for the molecule model specified in the URL and
%   displays the returned model. The URL must begin with the protocol, such
%   as http://, ftp:// and file://.
%
%   MOLVIEWER(PDBSTRUCT) loads a molecule model from a structure PDBSTRUCT.
%   PDBSTRUCT contains information of the molecule return by GETPDB or
%   PDBREAD function.
%
%   Examples:
%
%       % Display a model in a Elsevier MDL File.
%       h1 = molviewer('aspirin.mol');
%
%       % Display a model with PDB ID.
%       h2 = molviewer('2DHB');
%
%       % Load a H5N1 influenza virus hemagglutinin structure via URL.
%       h3 = molviewer('http://www.rcsb.org/pdb/files/2FK0.pdb.gz');
%
%       % Load a PDB structure
%       pdbstruct = getpdb('1vqx')
%       molviewer(pdbstruct)
%
%   Note: If you receive any errors related to memory or java heap space,
%   try increasing your java heap space as described here:
%
%         http://www.mathworks.com/support/solutions/data/1-18I2C.html
%
%   See also EVALRASMOLSCRIPT, GETPDB, MOLVIEWERDEMO, PDBDISTPLOT, PDBREAD,
%   PDBSUPERPOSE, PDBTRANSFORM, PDBWRITE, PROTEINPLOT, PROTEINPROPPLOT,
%   RAMACHANDRAN.

%   H = MOLVIEWER(FILENAME, HFIG) loads a molecule model file into the
%   viewer with handle HFIG.

% Copyright 2006-2012 The MathWorks, Inc.

% This function depends on a JVM.
if ~usejava('jvm')
    error(message('bioinfo:molviewer:NeedJVM', mfilename));
end

% start the viewer with no input
if nargin == 0
    hFig = initViewer([], true);
    resetViewer(hFig)
    if nargout > 0
        varargout{1} = hFig;
    end
    return;
end

callType = 'normal';

if nargin == 2
    hv = varargin{1};
    if isscalar(hv) && ishandle(hv) && strcmpi(get(hv, 'Tag'), 'BioinfoMolviewer')
        callType = 'hidden';
    else
        callType = 'normal';
        warning(message('bioinfo:molviewer:MoreThanOneInput'));
    end
    
elseif nargin > 2
    % Issue a warning and assume the call in 'normal'
    % which implies we use the first input argument and ignore all the
    % rest.
    warning(message('bioinfo:molviewer:MoreThanOneInput'));
end

% Check inputs
isPDBStruct = false;
ispdbid     = false;

switch callType
    case 'normal'
        switch class(mmsource)
            case 'struct'
                if isfield(mmsource, 'Model')
                    isPDBStruct = true;
                else
                    error(message('bioinfo:molviewer:InvalidInput'));
                end
            case 'char'
                mmsource = strtrim(mmsource);
                if size(mmsource, 1) == 1
                    if ~isempty(strfind(mmsource(1:min(10,end)), '://'))
                        % must be an url
                    elseif exist(mmsource, 'file') || exist(fullfile(pwd, mmsource), 'file')
                        % the file exists on MATLAB path
                        mmsource = handlefilename(mmsource);
                    else
                        % Check if we have a file with file extension
                        [pathstr,~,ext] = fileparts(mmsource);
                        if ~isempty(pathstr) || ~isempty(ext)
                            error(message('bioinfo:molviewer:InvalidPDBID',mmsource));
                        else
                            ispdbid = true;
                            mmsource = getPDBURL(mmsource);
                        end
                    end
                elseif isempty(mmsource)
                    %error(message('bioinfo:molviewer:InvalidInput'));
                else
                    error(message('bioinfo:molviewer:TooManyNames'));
                end
            otherwise
                error(message('bioinfo:molviewer:InvalidInput'));
        end
        
        hFig = initViewer([], false);
        if nargout > 0
            varargout{1} = hFig;
        end
        
        % For hidden option
    case 'hidden'
        hv = varargin{1};
        if isscalar(hv) && ishandle(hv) && strcmpi(get(hv, 'Tag'), 'BioinfoMolviewer')
            hFig = initViewer(hv, false);
            resetFigureTools(hFig);
            if strcmp(mmsource, 'resetFigureTool');
                if nargout > 0
                    varargout{1} = hFig;
                end
                return;
            end
        end
end

appdata = lGetAppData(hFig);
appdata.ispdbid = ispdbid;
appdata.mmsource = mmsource;

lSetAppData(hFig, appdata);

% Open the model
hFig = initFigureSize(hFig);

if isPDBStruct
    try       
        pstr = pdbwrite(mmsource);
    catch theException
        error(message('bioinfo:molviewer:FailReadPDBSTRUCT', theException.message));
    end
    loadModelFromString(hFig, pstr);
else
    loadModel(hFig, appdata.mmsource);
end

% Output handle
if nargout > 0
    varargout{1} = hFig;
end

end

%---------------- Callbacks ------------------------%
function doHelp(~, ~)
helpview(fullfile(docroot,'toolbox','bioinfo', 'bioinfo.map'), 'molviewer_refpage');
end

function doZoom(hSrc, hEvt,direction) %#ok
% Use the move command
move_script = sprintf('move 0 0 0 %s30 0 0 0 0 1', direction);
evalrasmolscript(gcbf, move_script);
end

function doSpin(hSrc, hEvt) %#ok
hfig=gcbf;
state = getToggleState(hfig, hSrc);

switch state
    case 'on'
        move_script = 'spin on';
    case 'off'
        move_script = 'spin off';
end
evalrasmolscript(hfig, move_script);
end

function doOpenConsole(hsrc, hevt) %#ok
appdata = lGetAppData(gcbf);
awtinvoke(appdata.viewer, 'openConsole()');
end

function doResetHomePosition(hSrc, hEvt)%#ok
appdata = lGetAppData(gcbf);
set(findall(appdata.toolbarbuttons,'Tag','spin'),'State','off')
set(findall(appdata.toolmitems,'Tag','spin_menu'),'Checked','off')
awtinvoke(appdata.viewer, 'resetHomePosition()');
end

function doBackground(hSrc, hEvt) %#ok
hfig=gcbf;
state = getToggleState(hfig, hSrc);

switch state
    case 'on'
        bg_script = 'background white';
    case 'off'
        bg_script = 'background black';
end
evalrasmolscript(hfig, bg_script);
end

function doSelect(hSrc, hEvt) %#ok
type = get(hSrc, 'Tag');
script = 'select none';
if strcmpi(type, 'select_all_menu')
    script = 'select all';
elseif strcmpi(type, 'select_none_menu')
    script = 'select none';
elseif strcmpi(type, 'select_protein_all_menu')
    script = 'select protein';
elseif strcmpi(type, 'select_protein_backbone_menu')
    script = 'select protein and backbone';
elseif strcmpi(type, 'select_protein_sidechain_menu')
    script = 'select protein and not backbone';
elseif strcmpi(type, 'select_nucleic_all_menu')
    script = 'select nucleic';
elseif strcmpi(type, 'select_nucleic_backbone_menu')
    script = 'select nucleic and backbone';
elseif strcmpi(type, 'select_nucleic_base_menu')
    script = 'select nucleic and not backbone';
elseif strcmpi(type, 'select_hetero_all_menu')
    script = 'select hetero';
elseif strcmpi(type, 'select_hetero_solvent_menu')
    script = 'select solvent';
elseif strcmpi(type, 'select_hetero_water_menu')
    script = 'select water';
elseif strcmpi(type, 'select_hetero_lifgand_menu')
    script = 'select ligand';
end
evalrasmolscript(gcbf, script)
end

function doLabels(hsrc, hEvt)%#ok
type = get(hsrc, 'Tag');
script = '';
appdata = lGetAppData(gcbf);
checkstates = {'on', 'off', 'off', 'off'};
if strcmpi(type, 'label_none_menu')
    script = 'label off';
elseif strcmpi(type, 'label_symbol_menu')
    if strcmp(get(hsrc, 'Checked'),'off')
        checkstates{2} = 'on';
        checkstates{1} = 'off';
        script = 'label %e';
    end
elseif strcmpi(type, 'label_name_menu')
    if strcmp(get(hsrc, 'Checked'),'off')
        checkstates{3} = 'on';
        checkstates{1} = 'off';
        script = 'label %a';
    end
elseif strcmpi(type, 'label_number_menu')
    if strcmp(get(hsrc, 'Checked'),'off')
        checkstates{4} = 'on';
        checkstates{1} = 'off';
        script = 'label %i';
    end
end

for i =1:4
    set(appdata.labelmitems(i), 'checked', checkstates{i});
end

if ~isempty(script)
    evalrasmolscript(gcbf, script)
end
end

function doShows(hsrc, hEvt)%#ok
type = get(hsrc, 'Tag');
chboxnum = [];

if strcmpi(type, 'axes_menu')
    chboxnum = 0;
elseif strcmpi(type, 'boundbox_menu')
    chboxnum = 1;
elseif strcmpi(type, 'unitcell_menu')
    chboxnum = 2;
elseif strcmpi(type, 'halos_menu')
    chboxnum = 3;
elseif strcmpi(type, 'dotsurface_menu')
    chboxnum = 4;
end
ison = isCheckOn(hsrc);
appdata = lGetAppData(gcbf);
awtinvoke(appdata.viewer, 'setShowCheckBox(IZ)', chboxnum, ison);
end

function ison = isCheckOn(hsrc)
ison = false;
if strcmp(get(hsrc, 'Checked'),'on')
    set(hsrc, 'Checked', 'off');
else
    set(hsrc, 'Checked', 'on');
    ison = true;
end
end

function showCallback(hsrc, hevt, hfig) %#ok
appdata = lGetAppData(hfig);
n = hevt.JavaEvent.item;
state = hevt.JavaEvent.state;
if ~isempty(appdata.showmitems)
    set(appdata.showmitems(n), 'Checked', char(state))
end
end

function notifyFileCallback(hsrc, hevt, hfig) %#ok
appdata = lGetAppData(hfig);
appdata.fileloaded = hevt.JavaEvent;
if appdata.fileloaded && appdata.noinputarg
    resetFigureTools(hfig);
    appdata.noinputarg = false;
end
lSetAppData(hfig, appdata)
drawnow
end


function modelInfoCallback(hsrc, hevt, hfig) %#ok
appdata = lGetAppData(hfig);
name = hevt.JavaEvent.name;
appdata.modelname = char(name);
if ~hevt.JavaEvent.isunitcell
    set(appdata.showmitems(3), 'Enable', 'off')
end

if ~hevt.JavaEvent.ispdb
    set(appdata.selectitems(1:3), 'Enable', 'off')
end
set(hfig, 'Name', [appdata.title ': ' appdata.modelname])
lSetAppData(hfig, appdata)
drawnow
end


function notifyFailedCallback(hsrc, hevt, hfig) %#ok
appdata = lGetAppData(hfig);

appdata.fileloaded = hevt.JavaEvent.loaded;
set(hfig, 'Name', [appdata.title,'Error Loading Model'])
lSetAppData(hfig, appdata)

drawnow
end

function doCopyImage(hsrc, hevt)%#ok
appdata = lGetAppData(gcbf);
awtinvoke(appdata.viewer, 'copyImageToClipboard()');
end

function doExport(hsrc, hevt)%#ok
type = get(hsrc, 'Tag');
imgtypes = {'JPEG', 'PNG', 'PPM'};
if strcmpi(type, 'exp_image_menu') ||strcmpi(type, 'Standard.SaveFigure')
    [filename, pathname, filteridx] = uiputfile(...
        {'*.jpg', [imgtypes{1} '(*.jpg)'];...
        '*.png', [imgtypes{2} '(*.png)'];...
        '*.ppm', [imgtypes{3} '(*.ppm)']},...
        'Save Image as');
    
    if ~filename
        return;
    end
    
    appdata = lGetAppData(gcbf);
    awtinvoke(appdata.viewer, 'saveViewImage(Ljava.lang.String;Ljava.lang.String;)',...
        fullfile(pathname, filename), imgtypes{filteridx});
elseif strcmpi(type, 'exp_pdf_menu')
    [filename, pathname] = uiputfile(...
        {'*.pdf', 'PDF(*.pdf)'},...
        'Save Image as', 'Untitled.pdf');
    
    if ~filename
        return;
    end
    appdata = lGetAppData(gcbf);
    awtinvoke(appdata.viewer, 'saveViewPDF(Ljava.lang.String;)',...
        fullfile(pathname, filename));
else
    return;
end
end

function doPrint(hsrc, hevt)%#ok
appdata = lGetAppData(gcbf);
awtinvoke(appdata.viewer, 'printView()');
end

function onViewerClosing(hfig, event)%#ok
appdata = lGetAppData(hfig);
if ~isempty(appdata.viewer)
    awtinvoke(appdata.viewer, 'cleanup()');
    delete(appdata.viewcontainer);
    delete(appdata.controlscontainer);
    delete(appdata.showlistener);
    delete(appdata.notifyfilelistener);
    delete(appdata.notifyfailedlistener);
    delete(appdata.modelinfolistener);
end
end

function resizeViewer(hfig, event) %#ok
set(hfig, 'Units', 'pixels');
resetViewer(hfig)
end

function resetViewer(hfig)
appdata = lGetAppData(hfig);
ppos = getpixelposition(hfig);
if ~isempty(appdata.viewcontainer)
    set(appdata.viewcontainer, 'position',[0 0 ppos(3) ppos(4)]);
end
drawnow
end

function figpos = updateViewerPosition(viewer, controller, hfig)
% Update viewer position according to the preferred size of the viewer
vd = awtinvoke(viewer, 'getPreferredSize()');
cd = awtinvoke(controller, 'getPreferredSize()');
width = vd.getWidth;
if cd.getWidth > width
    width = cd.getWidth;
end
ppos = getpixelposition(hfig);
ppos(3) = width;
ppos(4) = vd.getHeight;
setpixelposition(hfig,ppos);
figpos = get(hfig, 'Position');
end

%-------------Update toolbar and menu items --------------------------%
function appdata = createViewer(hfig)
appdata = lGetAppData(hfig);
try
    appdata.viewer = com.mathworks.toolbox.bioinfo.jmol.MolViewerPanel.getMolViewerPanel;
    
    if isempty(appdata.viewer) || appdata.viewer == 2
        error(message('bioinfo:molviewer:FailStartMolViewerPanel'));
    elseif appdata.viewer == 1
        error(message('bioinfo:molviewer:NoLicense'));
    end
    
    appdata.controller = awtinvoke(appdata.viewer, 'getControlPanel()');
    
    ppos = updateViewerPosition(appdata.viewer, appdata.controller, hfig);
    [appdata.viewerhandle, appdata.viewcontainer] =...
        javacomponent(appdata.viewer, [0 0 ppos(3) ppos(4)], hfig);
    
    % Add viewer mousewheellistener to Figure window.
    jf = getJavaFrame(hfig);
    wheelListeners = awtinvoke(appdata.viewer, 'getViewerMouseWheelListeners()');
    for i=1:wheelListeners.length
        jf.addMouseWheelListener(wheelListeners(i));
    end
    
    [appdata.controllerhandle, appdata.controlscontriner]= ...
        javacomponent(appdata.controller, java.awt.BorderLayout.SOUTH, handle(hfig));
    
    showcallback         = handle(appdata.viewer.getShowCallback);
    notifyfilecallback   = handle(appdata.viewer.getNotifyFileLoadedCallback());
    notifyfailedcallback = handle(appdata.viewer.getNotifyFileFailLoadedCallback());
    modelinfocallback    = handle(appdata.viewer.getModelInfoCallback());
    appdata.showlistener = handle.listener(showcallback, 'delayed', {@showCallback, hfig});
    appdata.notifyfilelistener = handle.listener(notifyfilecallback, 'delayed', {@notifyFileCallback, hfig});
    appdata.notifyfailedlistener = handle.listener(notifyfailedcallback,...
        'delayed', {@notifyFailedCallback, hfig});
    appdata.modelinfolistener = handle.listener(modelinfocallback,...
        'delayed', {@modelInfoCallback, hfig});
    
    appdata.isblankflag = false;
catch allExceptions %#ok<NASGU>
    onViewerClosing(hfig, []);
    delete(hfig)
    
    if isempty(appdata.viewer) || appdata.viewer == 2
        error(message('bioinfo:molviewer:FailInitializeMolViewer'));
    elseif appdata.viewer == 1
        error(message('bioinfo:molviewer:NoLicense'));
    else
        error(message('bioinfo:molviewer:ErrorOpeningViewer'));
    end
end
end

function initFigureTools(fig, isBlankFlag, varargin)
% helper function to set figure menus and toolbar
oldSH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

appdata = lGetAppData(fig);

appdata.isblankflag = isBlankFlag;

% Handle toolbar
if isempty(appdata.toolbar)
    appdata.toolbar = uitoolbar('parent', fig);
end
iconroot = fullfile(toolboxdir('bioinfo'), 'proteins', 'icons');
% Open icon should be on both cases
hOpen = uitoolfactory(appdata.toolbar, 'Standard.FileOpen');

iconfile = fullfile(iconroot, 'open_file.gif');
icon_cdata = lCreateToolbarIcon(iconfile);
set(hOpen, 'cdata', icon_cdata, 'ClickedCallback', @browseFile);

hSave = uitoolfactory(appdata.toolbar, 'Standard.SaveFigure');
iconfile = fullfile(iconroot, 'save_file.gif');
icon_cdata = lCreateToolbarIcon(iconfile);
set(hSave, 'cdata', icon_cdata,...
    'ClickedCallback', @doExport,...
    'tooltip', 'Export View Image');

hPrint = uitoolfactory(appdata.toolbar, 'Standard.PrintFigure');
iconfile = fullfile(iconroot, 'print.gif');
icon_cdata = lCreateToolbarIcon(iconfile);
set(hPrint,'cdata', icon_cdata, 'clickedcallback', @doPrint);

iconfile = fullfile(iconroot, 'zoom_in.gif');
hZoomin = createTool(appdata.toolbar, 1, iconfile,...
    'Zoom In', 'zoomin', 'separator', 'on');

iconfile = fullfile(iconroot, 'zoom_out.gif');
hZoomout = createTool(appdata.toolbar, 1, iconfile,...
    'Zoom Out', 'zoomout');

iconfile= fullfile(iconroot, 'spin_model.gif');
hSpin = createTool(appdata.toolbar, 2, iconfile,...
    'Spin On/Off',  'spin', 'State', 'off');

iconfile = fullfile(iconroot, 'background_black.gif');
hBg = createTool(appdata.toolbar,2, iconfile,...
    'Background B/W', 'background');

iconfile= fullfile(iconroot, 'reset_view.gif');
hReset = createTool(appdata.toolbar,1, iconfile,...
    'Reset Molecule to Home Position ','resetposition');

iconfile = fullfile(iconroot, 'scripting_console.gif');
hConsole = createTool(appdata.toolbar, 1, iconfile,...
    'Open Script Editor', 'console', 'separator', 'on');

set(hZoomin, 'clickedcallback', {@doZoom, ''});
set(hZoomout, 'clickedcallback', {@doZoom, '-'});
set(hReset, 'clickedcallback', {@doResetHomePosition});
set(hSpin, 'clickedcallback', {@doSpin});
set(hConsole, 'clickedcallback', {@doOpenConsole});
set(hBg, 'clickedcallback', {@doBackground});

iconfile = fullfile(iconroot, 'help.gif');
hHelp = createTool(appdata.toolbar, 1, iconfile,...
    'Help', 'help', 'separator', 'on');
set(hHelp, 'clickedcallback', {@doHelp});

appdata.toolbarbuttons = [hSave, hPrint, hZoomin,hZoomout, hSpin, hBg, hConsole, hReset];
% Hide the toolbar icons that are not needed for
if appdata.isblankflag
    for i = 1:length(appdata.toolbarbuttons)
        set(appdata.toolbarbuttons(i), 'visible', 'off');
    end
    set(appdata.toolbarbuttons(3), 'separator', 'off');
    set(appdata.toolbarbuttons(7), 'separator', 'off');
end

% Handle Menue items
% Delete figure menus not used
%h1 = findall(fig,'Type','uimenu', 'Label','&View');
h1 = findall(fig,'Type','uimenu', 'Tag','figMenuView');
%h2 = findall(fig,'Type','uimenu', 'Label','&Insert');
h2 = findall(fig,'Type','uimenu', 'Tag','figMenuInsert');
%h3 = findall(fig,'Type','uimenu', 'Label','&Desktop');
h3 = findall(fig,'Type','uimenu', 'Tag','figMenuDesktop');
delete([h1,h2, h3])

% Repair "File" menu
%hFileMenu = findall(fig,'Type','uimenu', 'Label','&File');
hFileMenu = findall(fig,'Type','uimenu', 'Tag','figMenuFile');
hf = get(hFileMenu,'children');

%h1 = findall(hFileMenu, 'Label', '&Open...');
h1 = findall(hFileMenu, 'Tag', 'figMenuOpen');
%h4 = findall(hFileMenu, 'Label', '&Close');
h4 = findall(hFileMenu, 'Tag', 'figMenuFileClose');
%h7 = findall(hFileMenu,'Label','&Print...');
h7 = findall(hFileMenu,'Tag','printMenu');
delete(setxor(hf,[h1,h4,h7]))

set(h1, 'Callback', @browseFile);
set(h4, 'Label', '&Close Viewer');
h2 = uimenu(hFileMenu, 'Label','Load PDB ID...', 'Position', 2, 'callback', @loadPDBID);%#ok
h3 = uimenu(hFileMenu,'Label','Open URL...', 'Position',3, 'callback', @loadURL);%#ok

h5 = createMenuItem(fig, hFileMenu, 'Export Image...', 'exp_image',...
    'separator', 'on',...
    'Position', 5,...
    'callback', @doExport);
h6 = createMenuItem(fig, hFileMenu, 'Export PDF...', 'exp_pdf',...
    'Position', 6, 'callback', @doExport);
set(h7, 'Callback', @doPrint);

appdata.filemenuitems = [h5,h6, h7];
if appdata.isblankflag
    for i = 1:length(appdata.filemenuitems)
        set(appdata.filemenuitems(i), 'visible', 'off')
    end
end

% Repair Edit menu
%hEditMenu = findall(fig,'Type','uimenu', 'Label','&Edit');
hEditMenu = findall(fig,'Type','uimenu', 'Tag','figMenuEdit');
hf = get(hEditMenu,'children');

%h1 = findall(hEditMenu, 'Label', 'Copy &Figure');
h1 = findall(hEditMenu, 'Tag', 'figMenuEditCopyFigure');
delete(setxor(hf,h1))
set(h1, 'callback', @doCopyImage);

createMenuItem(fig, hEditMenu, 'Select All', 'select_all',...
    'Position', 1,...
    'callback', @doSelect);
createMenuItem(fig, hEditMenu, 'Select None', 'select_none',...
    'Position', 2,...
    'callback', @doSelect);

% Selections
hsel1 = createMenuItem(fig, hEditMenu, 'Select Protein', 'select_protein',...
    'Position', 3,...
    'separator', 'on');
createMenuItem(fig, hsel1, 'All', 'select_protein_all',...
    'callback', @doSelect);
createMenuItem(fig, hsel1, 'Backbone', 'select_protein_backbone',...
    'callback',@doSelect);
createMenuItem(fig, hsel1, 'Side Chains', 'select_protein_sidechain',...
    'callback',@doSelect);

hsel2 = createMenuItem(fig, hEditMenu, 'Select Nucleic', 'select_nucleic',...
    'Position', 4);
createMenuItem(fig, hsel2, 'All', 'select_nucleic_all',...
    'callback', @doSelect);
createMenuItem(fig, hsel2, 'Backbone', 'select_nucleic_backbone',...
    'callback',@doSelect);
createMenuItem(fig, hsel2, 'Bases', 'select_nucleic_base',...
    'callback',@doSelect);

hsel3 = createMenuItem(fig, hEditMenu,'Select hetero', 'select_hetero',...
    'Position', 5);
createMenuItem(fig, hsel3, 'All', 'select_hetero_all',...
    'callback', @doSelect);
createMenuItem(fig, hsel3, 'Solvent', 'select_hetero_solvent',...
    'callback',@doSelect);
createMenuItem(fig, hsel3, 'Water', 'select_hectero_water',...
    'callback',@doSelect);
createMenuItem(fig, hsel3, 'Ligand', 'select_hectero_ligand',...
    'separator', 'on',...
    'callback',@doSelect);

hsel4 = createMenuItem(fig, hEditMenu,'Select Elements', 'select_element',...
    'Position', 6);
createMenuItem(fig, hsel4, 'All', 'select_hetero_all',...
    'callback', @doSelect);
createMenuItem(fig, hsel4, 'Solvent', 'select_hetero_solvent',...
    'callback',@doSelect);
createMenuItem(fig, hsel4, 'Water', 'select_hectero_water',...
    'callback',@doSelect);
appdata.selectitems = [hsel1, hsel2, hsel3, hsel4];

% Add Display menu
hDisplayMenu = uimenu(fig, 'Label','&Display',...
    'Tag','displaymenu',...
    'position', 3);

% Labels
labelitem = createMenuItem(fig, hDisplayMenu, 'Labels', 'labels',...
    'separator', 'on');
h1 = createMenuItem(fig, labelitem, 'None', 'label_none',...
    'checked', 'on',...
    'callback', @doLabels);
h2 = createMenuItem(fig, labelitem, 'Symbol', 'label_symbol',...
    'callback', @doLabels);
h3 = createMenuItem(fig, labelitem, 'Name', 'label_name',...
    'callback', @doLabels);
h4 = createMenuItem(fig, labelitem, 'Number', 'label_number',...
    'callback', @doLabels);

appdata.labelmitems =[h1, h2, h3, h4];
% Axes
hshow1 = createMenuItem(fig, hDisplayMenu, 'Axes', 'axes',...
    'separator', 'on',...
    'callback',  @doShows);
hshow2 = createMenuItem(fig, hDisplayMenu, 'Boundbox', 'boundbox',...
    'callback', @doShows);
hshow3 = createMenuItem(fig, hDisplayMenu, 'Unitcell', 'unitcell',...
    'callback', @doShows);
% Selection Halos
hshow4 = createMenuItem(fig, hDisplayMenu, 'Selection Halos', 'halos',...
    'separator', 'on',...
    'callback', @doShows);
% Dot surface
hshow5 = createMenuItem(fig, hDisplayMenu, 'Dot Surface', 'dotsurface',...
    'callback', @doShows);
appdata.showmitems = [hshow1,hshow2, hshow3, hshow4, hshow5];


% Repair Tools Menu
%hToolsMenu = findall(fig,'Type','uimenu', 'Label','&Tools');
hToolsMenu = findall(fig,'Type','uimenu','Tag','figMenuTools');
hf = get(hToolsMenu,'children');

delete(hf)
createMenuItem(fig, hToolsMenu, 'Zoom In', 'zoomin');
createMenuItem(fig, hToolsMenu, 'Zoom Out', 'zoomout');
h1 = createMenuItem(fig, hToolsMenu, 'Spin', 'spin');
h2 = createMenuItem(fig, hToolsMenu, 'Background', 'background',...
    'separator', 'on');
createMenuItem(fig, hToolsMenu, 'Reset Molecule Position', 'resetposition');
createMenuItem(fig, hToolsMenu, 'Scripting Console', 'console',...
    'separator', 'on');

appdata.toolmitems = [h1, h2];

% Repair Help menu
%hHelpMenu = findall(fig,'Type','uimenu','Label','&Help');
hHelpMenu = findall(fig,'Type','uimenu','Tag','figMenuHelp');
delete(get(hHelpMenu,'children'));
uimenu(hHelpMenu,'Label','Bioinformatics Toolbox Help',...
    'Position',1,...
    'Callback','helpview(fullfile(docroot,''toolbox'',''bioinfo'',''bioinfo.map''),''bioinfo_product_page'')')
uimenu(hHelpMenu,'Label','Molecule Viewer Help',...
    'Position',2,...
    'Callback',@doHelp)
uimenu(hHelpMenu,'Label','Examples',...
    'Position',3,...
    'Separator','on',...
    'Callback','demo(''toolbox'',''bioinfo'')')
tlbx = ver('bioinfo');
mailstr = ['web(''mailto:bioinfo-feedback@mathworks.com?subject=',...
    'Feedback%20for%20MolViewer%20in%20Bioinformatics',...
    '%20Toolbox%20',tlbx(1).Version,''')'];
uimenu(hHelpMenu,'Label','Send Feedback',...
    'Position',4,...
    'Separator','on',...
    'Callback',mailstr);

appdata.menus = [hEditMenu, hDisplayMenu, hToolsMenu];
if appdata.isblankflag
    for i = 1:length(appdata.menus)
        set(appdata.menus(i), 'visible', 'off');
    end
end

lSetAppData(fig, appdata);
set(0,'ShowHiddenHandles',oldSH)
end

function resetFigureTools(fig)
appdata = lGetAppData(fig);
for i = 1:length(appdata.toolbarbuttons)
    set(appdata.toolbarbuttons(i), 'visible', 'on');
end
set(appdata.toolbarbuttons(3), 'separator', 'on');
for i = 1:length(appdata.filemenuitems)
    set(appdata.filemenuitems(i), 'visible', 'on');
end

for i = 1:length(appdata.menus)
    set(appdata.menus(i), 'visible', 'on');
end

lSetAppData(fig, appdata);
end

function browseFile(hsrc, event)%#ok
[filename, pathname] = uigetfile(...
    {'*.pdb;*.ent;*.cif;*.mol;*.sdf;*.xyz;*.smol',...
    'Files(*.pdb, *.ent, *.cif, *.mol, *.sdf, *.xyz, *.smol)';...
    '*.*', 'All Files(*.*)'},...
    'Open');

if ~filename
    return;
else
    hfig = gcbf;
    mmsource = [pathname, filename];
    appdata = lGetAppData(hfig);
    if appdata.isblankflag
        molviewer(mmsource, hfig)
    else
        loadModel(hfig, mmsource);
    end
end
end

function loadPDBID(hsrc, evt) %#ok
prompt = {'Enter a PDB ID:'};
name = 'Input PDB ID';

answer = inputdlg(prompt, name);

if isempty(answer)
    return;
end

loadModel(gcbf, getPDBURL(answer));
end

function loadURL(hsrc, evt)%#ok
prompt = {'Enter URL of molecular model'};
name = 'Open URL';

answer = inputdlg(prompt, name);

if isempty(answer)
    return;
end
mmurl = answer{:};
if isempty(strfind(mmurl(1:min(10,end)), '://'))
    mmurl = ['http://', mmurl];
end
loadModel(gcbf, mmurl);
end

function loadModel(fig, mmsource)
appdata = lGetAppData(fig);
if isempty(appdata.viewer)
    molviewer(mmsource, fig);
else
    appdata.mmsource = mmsource;
    lSetAppData(fig, appdata)
    % Adding this switch for g524699
    if all(double(appdata.mmsource) < 128)
        awtinvoke(appdata.viewer, 'loadModelFile(Ljava.lang.String;)', appdata.mmsource);
    else
        awtinvoke(appdata.viewer, 'evalMVScriptT(Ljava/lang/String;)', sprintf('load %s', appdata.mmsource));
    end
    drawnow;
end
end

% If we have a model file in memory, this function passes the model to Jmol
% directly.
function loadModelFromString(hFig, modelString)
appdata = lGetAppData(hFig);
awtinvoke(appdata.viewer, 'loadInline(Ljava/lang/String;)', modelString);
end

function pdburl = getPDBURL(pdbid)
rcsbPDBURLHeader = 'http://www.rcsb.org/pdb/files/';
rcsbPDBURLTail = '.pdb.gz';

pdburl = [rcsbPDBURLHeader char(pdbid) rcsbPDBURLTail];
end

function hfig = initFigureSize(hfig)
% Resize the whole figure window
appdata = lGetAppData(hfig);
units = get(hfig, 'Units');
set(hfig, 'Units', 'pixels');

pos = getpixelposition(hfig);
pos2 = [pos(1) pos(2) pos(3)  pos(4)];
if ismac
	pos2(3) = 680;
    pos2(4) = 460;
end
drawnow
set(hfig, 'position', pos2);
set(hfig, 'units', units')

lSetAppData(hfig, appdata);
resizeViewer(hfig,[]);
end

function hfig = initViewer(hfig, isBlankFlag)
if isempty(hfig)
    hfig = figure('WindowStyle', 'normal', ...
        'Tag', 'BioinfoMolviewer',...
        'toolbar', 'none',...
        'units', 'pixels',...
        'DeleteFcn', @onViewerClosing,...
        'resize', 'on',...
        'ResizeFcn', @resizeViewer,...
        'PaperPositionMode', 'auto', ...
        'Name', 'Molecule Viewer: Blank',...
        'numbertitle','off', ...
        'visible', 'on');
    initFigureTools(hfig, isBlankFlag)
end

appdata = createViewer(hfig);

if isBlankFlag
    appdata.noinputarg = true;
end
lSetAppData(hfig, appdata);

end

function lSetAppData(hfig,appdata)
setappdata(hfig,'BioMolViewer',appdata);
end

function appdata = lGetAppData(hfig)
if isappdata(hfig,'BioMolViewer')
    appdata = getappdata(hfig,'BioMolViewer');
else
    appdata = guihandles(hfig);
    appdata.ispdbid = true;
    appdata.mmsource = [];  % molecule model source
    appdata.viewer = [];
    appdata.viewerhandle = [];
    appdata.viewcontainer = [];
    appdata.controller = [];
    appdata.controllerhandle=[];
    appdata.controlscontainer = [];
    appdata.loadstatus = '';
    
    appdata.isblankflag = false;
    appdata.toolbar = [];
    appdata.toolbarbuttons = [];
    appdata.filemenuitems = [];
    appdata.showmitems = [];
    appdata.selectitems = [];
    appdata.labelmitems = [];
    appdata.toolmitems = [];
    appdata.menus = [];
    
    appdata.showlistener = [];
    appdata.modelinfolistener = [];
    appdata.notifyfilelistener = [];
    appdata.notifyfailedlistener = [];
    appdata.title = 'Molecule Viewer';
    appdata.modelname = [];
    appdata.fileloaded = false;
    appdata.noinputarg = false;
end
end
%-------- Helper functions-----------------------------------
function state = getToggleState(hfig, hSrc)
type = get(hSrc,'type');
appdata = lGetAppData(hfig);
if strcmpi(type,'uimenu')
    tag = strrep(get(hSrc,'tag'),'_menu','');
    hSrcPeer = findall(hfig,'tag', tag);
    
    if strcmp(get(hSrcPeer,'State'),'on')
        set(hSrcPeer,'State','off');
    else
        set(hSrcPeer,'State','on');
    end
    state = get(hSrcPeer, 'State');
    set(hSrc, 'Checked', state);
elseif strcmpi(type, 'uitoggletool')
    menu_tag = [get(hSrc,'tag'),'_menu'];
    hSrcPeer = findobj(appdata.toolmitems, 'tag',menu_tag);
    state = get(hSrc,'State');
    
    if ~isempty(hSrcPeer)
        set(hSrcPeer,'Checked', state);
    end
end
end

function menuH = createMenuItem(fig, parent, labelname, tag, varargin)
% CREATEMENUITEM Creates a menu Item and connects its behavior
%    to the behavior of other ui components with the same tag.

% Ceate File Menu etc.
if isempty(parent)
    menuH = uimenu(fig,'Label',labelname,'Tag', tag, varargin{:});
    return;
end

menu_tag = [tag, '_menu'];
menuH = uimenu(parent,'Label',labelname,'Tag', menu_tag, varargin{:});

comp = findobj(fig,'tag',tag);
if ~isempty(comp)
    if strcmpi(get(comp,'type'),'uipushtool') % for toolbar pushbuttons
        set(menuH,'Callback',get(comp,'ClickedCallback'));
    elseif strcmpi(get(comp,'type'),'uitoggletool') % for toolbar toggles
        set(menuH,'Callback',get(comp,'ClickedCallback'));
    elseif strcmpi(get(comp,'type'),'uicontrol') %  for marker push buttons
        set(menuH,'Callback',get(comp,'Callback'));
    end
end
end

function htool = createTool(toolbar, type, iconfile, tooltip, tag, varargin)
% create toolbar buttons% iType = 1 - uipushtool, 2-uitoggletool

icon_cdata = lCreateToolbarIcon(iconfile);
props = {toolbar,'cdata', icon_cdata, 'Tooltip', tooltip, 'Tag', tag, varargin{:}}; %#ok<CCAT>
switch type
    case 1
        htool=uipushtool(props{:});
    case 2
        htool=uitoggletool(props{:});
end
end

function cdata = lCreateToolbarIcon(filename)
% LOCALCREATETOOLBARICON read an image file and convert it to CData for a
% toolbar icon.
% Adopted from Bill York's iconread.m

[~,~,ext] = fileparts(filename);
% if this is a mat-file, look for the variable cdata (or something like it)
if isequal(lower(ext),'.mat')
    data = load(filename,'cdata');
    cdata = data.cdata;
    return
end

[cdata,map] = imread(filename);
if isempty(cdata)
    return;
end

if isempty(map)
    % need to use doubles, nan's only work as doubles
    cdata = double(cdata);
    cdata = cdata/255;
else
    % Set all white (1,1,1) colors to be transparent (nan)
    ind = find(map(:,1)+map(:,2)+map(:,3)==3);
    map(ind) = nan; %#ok
    cdata = ind2rgb(cdata,map);
end
end

function fullname = handlefilename(filename)
wfile = which(filename);
if ~isempty(wfile)
    filename = wfile;
end

[thePath, theName, theExt] = fileparts(filename);
if isempty(thePath)
    thePath = pwd;
end

fullname = fullfile(thePath, [theName theExt]);

if strncmp(fullname, '\\', 2)
    % Convert network paths to file:// URIs to work around g1232432
    fullname = ['file://' strrep(fullname, '\', '/')];
end
end

function javaFrame = getJavaFrame(f)
% store the last warning thrown
[ lastWarnMsg lastWarnId ] = lastwarn;

% disable the warning when using the 'JavaFrame' property
% this is a temporary solution
oldJFWarning = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
javaFrame = get(f,'JavaFrame');
warning(oldJFWarning.state, 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

% restore the last warning thrown
lastwarn(lastWarnMsg, lastWarnId);
end


