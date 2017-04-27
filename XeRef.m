classdef XeRef < handle
    
    properties
        
        data
        layers
        
        protein
        pdbFiles
        
        gui
        handles
        
    end
    
    methods
        
        function this = XeRef(files)
            
            createView();
            createController();
            this.control('initialize');
            if nargin == 1
                this.control('load-data', files);
            end
            
            function createView()
                
                initializeView();
                createFilePanel();
                createAxes();
                createFigureControls();
                createRightPanel();
                
                createBasicInfoTalbe();
                createParameterTable();
                
                createFittingControls();
                createOutputandSaveButtons();
                
                function initializeView()
                    
                    screenSize = [335, 300, 1250, 800];
                    
                    set(0, 'units', 'pixels');
                    pix = get(0, 'screensize');
                    if pix(4) * 0.85 <= 800
                        screenSize(4) = pix(4)*0.85;
                    end
                    
                    this.handles = figure('Visible', 'off', 'Name', 'XeRef', 'NumberTitle', 'off', 'Units', 'pixels',...
                        'Position', screenSize, 'Resize', 'on');
                    movegui(this.handles, 'center');
                    this.handles.Visible = 'on';
                    
                end
                
                function createFilePanel()
                    
                    panel = uipanel(this.handles, 'Title', 'Data Files', 'Units', 'normalized',...
                        'Position', [0.014 0.02 0.16 0.97]);
                    
                    uicontrol(panel,'Style','text','String','Reflectivity Scans','Units','normalized',...
                        'Position',[0.05 0.965 0.8 0.03]);
                    
                    this.gui.dataFiles = uicontrol(panel,'Style','listbox','Units','normalized',...
                        'Position',[0.05 0.56 0.9 0.405],'Max',2);
                    
                    this.gui.loadData = uicontrol(panel,'Style','pushbutton','String','Load','Units','normalized',...
                        'Position',[0.035 0.52 0.3 0.032]);
                    
                    this.gui.deleteData = uicontrol(panel,'Style','pushbutton','String','Delete','Units','normalized',...
                        'Position',[0.38 0.52 0.3 0.032]);
                    
                    uicontrol(panel,'Style','text','String','Protein structures','Units','normalized',...
                        'Position',[0.05 0.49 0.8 0.03]);
                    
                    this.gui.pdbFiles = uicontrol(panel,'Style','listbox','Units','normalized',...
                        'Position',[0.05 0.051 0.9 0.44],'Max',2);
                    
                    this.gui.loadPdb = uicontrol(panel,'Style','pushbutton','String','Load PDB','Units','normalized',...
                        'Position',[0.035 0.015 0.35 0.032]);
                    
                    this.gui.deletePdb = uicontrol(panel,'Style','pushbutton','String','Delete','Units','normalized',...
                        'Position',[0.38 0.015 0.3 0.032]);
                    
                end
                
                function createAxes()
                    
                    ax1 = axes('Parent',this.handles,'Units','normalized','Position',[0.215 0.52 0.38 0.44]);
                    ax1.XLim = [0 10];
                    ax1.YLim = [0 10];
                    ax1.XTick = [0 2 4 6 8 10];
                    ax1.YTick = [0 2 4 6 8 10];
                    ax1.XLabel.String = 'x1';
                    ax1.YLabel.String = 'y1';
                    
                    this.gui.ax1 = ax1;
                    
                    % plot region 2
                    ax2 = axes('Parent',this.handles,'Units','normalized','Position',[0.215 0.08 0.38 0.35]);
                    ax2.XLim = [0 10];
                    ax2.YLim = [0 10];
                    ax2.XTick = [0 2 4 6 8 10];
                    ax2.YTick = [0 2 4 6 8 10];
                    ax2.XLabel.String = 'x2';
                    ax2.YLabel.String = 'y2';
                    
                    this.gui.ax2 = ax2;
                    
                end
                
                function createFigureControls()
                    
                    this.gui.normalized = uicontrol(this.handles, 'Style','checkbox','String','Normalized','Units','normalized','Visible','on',...
                        'Position',[0.6 0.94 0.1 0.018], 'Value', 1);
                    
                    this.gui.showCal = uicontrol(this.handles,'Style','checkbox','String','Show Calc.','Units','normalized',...
                        'Position',[0.6 0.915 0.1 0.018], 'Value', 1);
                    
                    this.gui.showFit = uicontrol(this.handles,'Style','checkbox','String','Show Fit','Units','normalized',...
                        'Position',[0.6 0.89 0.1 0.018], 'Value', 0);
                    
                    this.gui.showEd = uicontrol(this.handles,'Style','checkbox','String','Ed Profile','Units','normalized',...
                        'Position',[0.6 0.41 0.1 0.018], 'Value', 1);
                    
                    this.gui.showProtein = uicontrol(this.handles,'Style','checkbox','String','Protein','Units','normalized',...
                        'Position',[0.6 0.385 0.1 0.018], 'Value', 0);
                    
                end
            
                function createRightPanel()
                    
                    this.gui.rightPanel = uipanel(this.handles, 'Units','normalized','Position',[0.68 0.02 0.31 0.97]);
                    
                end
                
                function createBasicInfoTalbe()
                    
                    rightPanel = this.gui.rightPanel;
                    
                    rowName = {'Beam Energy (keV)'; 'Roughness (A)'; 'Cut Qz Below'};
                    colName = {};
                    columnFormat = {'numeric'};
                    columnWidth = {120};
                    tableData = {10; 3.4; 0.03};
                    
                    this.gui.basicInfoTable = uitable(rightPanel, 'Data', tableData, 'ColumnName', colName, ...
                        'ColumnFormat', columnFormat, 'ColumnEditable', true, 'Units','normalized', ...
                        'ColumnWidth',columnWidth,'RowName',rowName, 'RowStriping','off',...
                        'Position', [0.025 0.91 0.935 0.08], 'TooltipString', 'Press enter to update value.');
                    
                end
                
                function createParameterTable()
                    
                    rightPanel = this.gui.rightPanel;
                    
                    rowName = {'Qz-Offset', 'Top-ED', 'Layer-2-ED', 'Layer-1-ED', 'Bottom-Ed', 'Layer-2-Thkns', 'Layer-1-Thkns'};
                    colName = {'Min','Max','Start','Fix','Plot'};
                    colFormat = {'numeric','numeric','numeric','logical','logical'};
                    colWidth = {55 55 55 30 30};
%                     tableData = {-0.01, 0.01, 0, false, false; 0, 0, 0, true, false; 0.2, 0.3, 0.26, false, false;...
%                         0.4, 0.5, 0.44, false, false; 0.335, 0.335, 0.335, true, false; 8, 16, 12, false, false;...
%                         6, 14, 10, false, false};
                    tableData = {-5e-4, 5e-4, 0, false, false; 0, 0, 0, true, false; 0.18, 0.32, 0.24, false, false;...
                        0.38, 0.48, 0.46, false, false; 0.335, 0.335, 0.335, true, false; 9, 15, 12, false, false;...
                        8, 20, 12, false, false};
                    
                    h = 0.88;
                    
                    this.gui.parametersTableTitle = uicontrol(rightPanel,'Style','text','String','Fitting Parameters:',...
                        'Units','normalized','HorizontalAlignment','left', 'Position', [0.025 h 0.3 0.025]);
                    
                    this.gui.toggleProtein = uicontrol(rightPanel,'Style','radiobutton','String','Add Protein Layer',...
                        'Units','normalized', 'Position',[0.68 h 0.4 0.03]);
                    
                    this.gui.parametersTable = uitable(rightPanel,'Data', tableData,'ColumnName', colName,...
                        'ColumnFormat', colFormat,'ColumnEditable', [true true true true true],'Units','normalized',...
                        'ColumnWidth',colWidth,'RowName',rowName,'RowStriping','off',...
                        'Position', [0.025 0.455 0.935 0.42]);
                    
                    this.gui.addLayer = uicontrol(rightPanel,'Style','pushbutton','String', 'Add', 'Units','normalized',...
                        'Position', [0.725 0.42 0.11 0.03]);
                    
                    this.gui.deleteLayer = uicontrol(rightPanel,'Style','pushbutton','String', 'Delete',...
                        'Units','normalized', 'Position', [0.84 0.42 0.12 0.03]);
                    
                end
                
                function createFittingControls()
                    
                    rightPanel = this.gui.rightPanel;
                    
                    h = 0.42;
                    
                    this.gui.loadPara = uicontrol(rightPanel,'Style','pushbutton','String','Load Para','Units','normalized',...
                        'Position',[0.024 h-0.03 0.17 0.03]);
                    
                    this.gui.savePara = uicontrol(rightPanel,'Style','pushbutton','String','Save Para','Units','normalized',...
                        'Position',[0.19 h-0.03 0.17 0.03]);
                    
                    this.gui.stepInput = uicontrol(rightPanel,'Style','edit','String', 10, 'Units','normalized',...
                        'HorizontalAlignment','left','Position',[0.4 h-0.03 0.1 0.03]);
                    
                    this.gui.stepText = uicontrol(rightPanel,'Style','text','String','Steps','Units','normalized',...
                        'HorizontalAlignment','left','Position', [0.515 h-0.035 0.08 0.03]);
                    
                    this.gui.angleFit = uicontrol(rightPanel,'Style','pushbutton','String', 'Angle Fit',...
                        'Units','normalized', 'Position', [0.024 0.42 0.15 0.03]);
                    
                    this.gui.fitButton = uicontrol(rightPanel,'Style','pushbutton','String','Fit','Units','normalized',...
                        'Position',[0.82 h-0.03 0.15 0.03]);
                    
                    this.gui.quickFit = uicontrol(rightPanel,'Style','pushbutton','String', 'Quick Fit',...
                        'Units','normalized', 'Position', [0.65 h-0.03 0.15 0.03]);
                    
                    this.gui.withText = uicontrol(rightPanel,'Style','text','String','With','Units','normalized','HorizontalAlignment','left',...
                        'Position',[0.025 h-0.065 0.07 0.03]);
                    
                    this.gui.confidenceInput = uicontrol(rightPanel,'Style','edit','String','95','Units','normalized',...
                        'HorizontalAlignment','left','Position',[0.1 h-0.06 0.07 0.03]);
                    
                    this.gui.confidenceText = uicontrol(rightPanel,'Style','text','String','% confidence window','Units','normalized','HorizontalAlignment','left',...
                        'Position',[0.171 h-0.065 0.28 0.03]);
                    
                    this.gui.recordFitting = uicontrol(rightPanel,'Style','pushbutton','String','Record Fitting','Units','normalized',...
                        'Position',[0.452 h-0.06 0.22 0.03]);
                    
                    this.gui.updateStartButton = uicontrol(rightPanel,'Style','pushbutton','String','Update Starts','Units','normalized',...
                        'Position',[0.75 h-0.06 0.22 0.03]);
                    
                end
                
                function createOutputandSaveButtons()
                    
                    rightPanel = this.gui.rightPanel;
                    
                    this.gui.output = uicontrol(rightPanel,'Style','edit','Max',2,'HorizontalAlignment','left','Units','normalized',...
                        'Position',[0.03 0.07 0.935 0.28]);
                    
                    this.gui.clearOutput = uicontrol(rightPanel,'Style','pushbutton','String','Clear','Units','normalized',...
                        'Position',[0.82 0.038 0.15 0.03]);
                    
                    uicontrol(rightPanel,'Style','text','String','Save:','Units','normalized',...
                        'HorizontalAlignment','left','Position',[0.025 0.035 0.08 0.025]);
                    
                    this.gui.saveOutput = uicontrol(rightPanel,'Style','pushbutton','String','Output Text','Units','normalized',...
                        'Position',[0.024 0.007 0.2 0.03]);
                    
                    this.gui.saveUpperFigure = uicontrol(rightPanel,'Style','pushbutton','String','Upper Figure','Units','normalized',...
                        'Position',[0.234 0.007 0.2 0.03]);
                    
                    this.gui.saveLowerFigure = uicontrol(rightPanel,'Style','pushbutton','String','Lower Figure','Units','normalized',...
                        'Position',[0.444 0.007 0.2 0.03]);
                    
                    this.gui.saveData = uicontrol(rightPanel,'Style','pushbutton','String','Data & Fit','Units','normalized',...
                        'Position',[0.66 0.007 0.17 0.03]);
                    
                end
                
            end
            
            function createController()
                
                % files
                this.gui.loadData.Callback = @(varargin) this.control('load-data');
                this.gui.dataFiles.Callback = @(varargin) this.control('choose-data');
                this.gui.loadPdb.Callback = @(varargin) this.control('load-pdb');
                this.gui.pdbFiles.Callback = @(varargin) this.control('choose-pdb');
                
                % figure controls
                this.gui.normalized.Callback = @(varargin) this.control('toggle-fresnel');
                this.gui.showCal.Callback = @(varargin) this.control('toggle-cal');
                this.gui.showFit.Callback = @(varargin) this.control('toggle-fit');
                this.gui.showProtein.Callback = @(varargin) this.control('show-protein');
                this.gui.showEd.Callback = @(varargin) this.control('show-ed');
                
                % table controls
                this.gui.basicInfoTable.CellEditCallback = @(source, eventdata, varargin) this.control('basic-info', eventdata);
                this.gui.toggleProtein.Callback = @(varargin) this.control('toggle-protein');
                this.gui.addLayer.Callback = @(varargin) this.control('add-layer');
                this.gui.deleteLayer.Callback = @(varargin) this.control('delete-layers');
                this.gui.parametersTable.CellEditCallback = @(source, eventdata, varargin) this.control('parameter-table', eventdata);
                
                % fitting
                this.gui.fitButton.Callback = @(varargin) this.control('fit');
                this.gui.updateStartButton.Callback = @(varargin) this.control('update-starts');
                this.gui.quickFit.Callback = @(varargin) this.control('quick-fit');
                this.gui.angleFit.Callback = @(varargin) this.control('angle-fit');
                
            end
            
        end
        
        function control(this, trigger, varargin)
            
            state = guiState();
            
            switch state
                case 'empty'
                    switch trigger
                        case 'initialize'
                            pData = readParameterTable();
                            this.model(state, trigger, pData);
                            this.view(state, trigger);
                        case 'load-data'
                            this.model(state, trigger, varargin{1});
                            this.view(state, trigger);
                        otherwise
                            sprintf('State: %s, trigger: %s is not found for the controller', state, trigger)
                    end
                case 'inspect'
                    switch trigger
                        case 'load-data'
                            this.model(state, trigger);
                            this.view(state, trigger);
                        case 'load-pdb'
                            if ~isempty(varargin) > 0
                                pdbfiles = varargin{1};
                            else
                                pdbfiles = [];
                            end
                            this.model(state, trigger, pdbfiles);
                            this.view(state, trigger); 
                        case 'choose-pdb'
                            this.model(state, trigger);
                            this.view(state, trigger);
                        case {'choose-data', 'toggle-fresnel', 'toggle-cal'}
                            this.view(state, trigger);
                        case 'show-protein'
                            this.view(state, trigger);
                        case 'show-ed'
                            this.view(state, trigger);
                        case 'parameter-table'
                            eventdata = varargin{1};
                            [viewUpate, modelUpdate] = paraTableEditEntailsUpdate(eventdata);
                            paras = readParameterTable();
                            if modelUpdate
                                this.model(state, trigger, paras);
                            end
                            if viewUpate
                                this.view(state, trigger, paras.theta, paras.phi);
                            end
                        case 'basic-info'
                            eventdata = varargin{1};
                            needUpdate = handleBasicInfoEdit(eventdata);
                            if needUpdate
                                t = readBasicTableData();
                                this.model(state, trigger, t);
                                this.view(state, trigger, t);
                            end
                        case 'fit'
                            refData = this.data{this.gui.dataFiles.Value(1)};
                            pData = readParameterTable();
                            steps = str2double(this.gui.stepInput.String);
                            this.model(state, trigger, refData, pData, steps);
                            this.view(state, trigger);
                        case 'quick-fit'
                            refData = this.data{this.gui.dataFiles.Value(1)};
                            pData = readParameterTable();
                            steps = str2double(this.gui.stepInput.String);
                            this.model(state, trigger, refData, pData, steps);
                            this.view(state, trigger);
                        case 'angle-fit'
                            refData = this.data{this.gui.dataFiles.Value(1)};
                            pData = readParameterTable();
                            steps = str2double(this.gui.stepInput.String);
                            this.model(state, trigger, refData, pData, steps);
                            this.view(state, trigger);
                        case 'toggle-fit'
                            if isfield(this.layers.fits, 'all')
                                this.view(state, trigger);
                            end
                        case 'update-starts'
                            newStarts = this.layers.fits.all.para_all;
                            this.view(state, trigger, newStarts);
                            paras = readParameterTable();
                            this.model(state, trigger, paras);
                            this.view(state, 'update-starts-follow-up');
                        case 'toggle-protein'
                            if isempty(this.gui.pdbFiles.String)
                                this.gui.toggleProtein.Value = false;
                            else
                                this.view(state, trigger, this.gui.toggleProtein.Value);
                                paras = readParameterTable();
                                this.model(state, trigger, paras);
                                this.view(state, 'toggle-protein-update');
                            end
                        otherwise
                            sprintf('State: %s, trigger: %s is not found for the controller', state, trigger)
                    end
                otherwise
                    sprintf('State: %s not found for the controller', state)
            end
            
            function state = guiState()
                
                if isempty(this.gui.dataFiles.String)
                    state = 'empty';
                else
                    state = 'inspect';
                end
                
            end
            
            function [viewUpdate, modelUpdate] = paraTableEditEntailsUpdate(eventdata)
                
                viewUpdate = false;
                modelUpdate = false;
                ind1 = eventdata.Indices(1);
                ind2 = eventdata.Indices(2);
                newData = eventdata.NewData;
                table = this.gui.parametersTable;
                
                switch ind2
                    case 1
                        if newData > table.Data{ind1, 2}
                            table.Data{ind1, 2} = newData;
                            table.Data{ind1, 3} = newData;
                            table.Data{ind1, 4} = true;
                            viewUpdate = true;
                            modelUpdate = true;
                        elseif newData > table.Data{ind1, 3}
                            table.Data{ind1, 3} = newData;
                            viewUpdate = true;
                            modelUpdate = true;
                        else
                            table.Data{ind1, 4} = false;
                        end
                    case 2
                        if newData < table.Data{ind1, 1}
                            table.Data{ind1, 1} = newData;
                            table.Data{ind1, 3} = newData;
                            table.Data{ind1, 4} = true;
                            viewUpdate = true;
                            modelUpdate = true;
                        elseif newData < table.Data{ind1, 3}
                            table.Data{ind1, 3} = newData;
                            viewUpdate = true;
                            modelUpdate = true;
                        else
                            table.Data{ind1, 4} = false;
                        end
                    case 3
                        if newData < table.Data{ind1, 1}
                            table.Data{ind1, 1} = newData;
                            table.Data{ind1, 4} = false;
                        elseif newData > table.Data{ind1, 2}
                            table.Data{ind1, 2} = newData;
                            table.Data{ind1, 4} = false;
                        end
                        viewUpdate = true;
                        modelUpdate = true;
                    case 4
                        if table.Data{ind1, ind2}
                            table.Data{ind1, 1} = table.Data{ind1, 3};
                            table.Data{ind1, 2} = table.Data{ind1, 3};
                        else
                            switch ind1
                                case 1
                                    table.Data{ind1, 1} = -0.001;
                                    table.Data{ind1, 2} = 0.001;
                                case 2
                                    table.Data{ind1, ind2} = true;
                                    disp('Does not support a different top phase yet.');
                                otherwise
                                    table.Data{ind1, 1} = 0.5 * table.Data{ind1, 3};
                                    table.Data{ind1, 2} = 2 * table.Data{ind1, 3};
                            end
                        end
                    case 5
                        if isfield(this.layers.fits, 'one') && this.layers.fits.all.fitted(ind1)
                            sel = cell2mat(table.Data(:, end));
                            num = sum(sel);
                            if num >= 3
                                indices = find(sel);
                                indices = indices(indices ~= ind1);
                                table.Data{indices(1), ind2} = false;
                            end
                            viewUpdate = true;
                        else
                            table.Data{ind1, ind2} = false;
                        end
                end
                
            end
            
            function t = readParameterTable()
                
                dat = this.gui.parametersTable.Data;
                [m, ~] = size(dat);
                
                mat = cell2mat(dat(:, 1:3));
                
                t.p0 = mat(:, 3)';
                t.lb = mat(:, 1)';
                t.ub = mat(:, 2)';
                
                t.pro = this.gui.toggleProtein.Value;
                
                if t.pro
                    n_layer = (m - 4) / 2;
                    t.phi = t.p0(end - 1);
                    t.theta = t.p0(end - 2);
                    t.insertion = t.p0(end - 3);
                    t.density = t.p0(end - 4);
                    t.protein = this.protein;
                    t.inhomo = t.p0(end);
                else
                    n_layer = (m + 1) / 2;
                end
                
                t.ed = t.p0(n_layer + 1 : -1 : 2);
                t.thickness = [Inf, t.p0(n_layer * 2 - 1 : -1 : 2 + n_layer), Inf];
                
            end
            
            function needUpdate = handleBasicInfoEdit(eventdata)
                
                needUpdate = false;
                
                ind1 = eventdata.Indices(1);
                
                switch ind1
                    case 1
                    case 2
                        needUpdate = true;
                    case 3
                end
                
            end
            
            function t = readBasicTableData()
                
                t = cell2mat(this.gui.basicInfoTable.Data)';
                
            end
            
        end
        
        function model(this, state, trigger, varargin)
            
            switch state
                case 'empty'
                    switch trigger
                        case 'initialize'
                            pData = varargin{1};
                            this.layers = RefLayers(10, pData.ed, pData.thickness);
                        case 'load-data'
                            if nargin == 4
                                datafiles = varargin{1};
                                loadNewData(datafiles);
                            else
                                loadNewData();
                            end
                    end
                case 'inspect'
                    switch trigger
                        case 'load-data'
                            loadNewData();
                        case 'load-pdb'
                            pdbfiles = varargin{1};
                            loadNewPdb(pdbfiles);
                        case 'layer-table'
                            dat = varargin{1};
                            this.layers = RefLayers(10, dat.ed, dat.thickness);
                        case 'parameter-table'
                            paras = varargin{1};
                            this.layers.updateModel(paras);
                        case 'basic-info'
                            tdata = varargin{1};
                            this.layers.updateRoughness(tdata(2));
                        case 'fit'
                            refData = varargin{1};
                            pData = varargin{2};
                            steps = varargin{3};
                            this.layers.fitDataThorough(refData, pData.p0, pData.lb, pData.ub, steps);
                        case 'quick-fit'
                            refData = varargin{1};
                            pData = varargin{2};
                            steps = varargin{3};
                            this.layers.fitDataQuick(refData, pData.p0, pData.lb, pData.ub, steps);
                        case 'angle-fit'
                            refData = varargin{1};
                            pData = varargin{2};
                            steps = varargin{3};
                            this.layers.fitAngleGrid(refData, pData.p0, pData.lb, pData.ub, steps);
                        case 'update-starts'
                            paras = varargin{1};
                            this.layers.updateModel(paras);
                        case 'choose-pdb'
                            this.protein = RefProtein(this.pdbFiles{this.gui.pdbFiles.Value(1)});
                        case 'toggle-protein'
                            paras = varargin{1};
                            this.layers.updateModel(paras);
                    end
                otherwise
                    sprintf('Case: %s is not found for the view', state);
            end
            
            function loadNewData(datafiles)
                
                if nargin == 1
                    files = datafiles;
                    newdata = cell(1, length(files));
                    
                    for i = 1 : length(files)
                        newdata{i} = RefData(files{i});
                    end
                    
                    this.data = [this.data, newdata];
                else
                    [files, path] = uigetfile('*.ref', 'Select the data files.', 'MultiSelect', 'on');
                    if ~isnumeric(files)
                        if isa(files, 'char')
                            files = {files};
                        end
                        
                        newdata = cell(1, length(files));
                        
                        for i = 1 : length(files)
                            newdata{i} = RefData(fullfile(path, files{i}));
                        end
                        
                        this.data = [this.data, newdata];
                    end
                end
                
            end
            
            function loadNewPdb(files)
                
                if nargin == 0 || isempty(files)
                    [files, path] = uigetfile('*.pdb', 'Select the pdb files.', 'MultiSelect', 'on');
                    if ~isnumeric(files)
                        if isa(files, 'char')
                            files = {files};
                        end
                        for i = 1 : length(files)
                            files{i} = fullfile(path, files{i});
                        end
                    end
                elseif isa(files, 'char')
                    files = {files};
                end
                
                this.pdbFiles = [this.pdbFiles, files];
                
                if isempty(this.protein)
                    this.protein = RefProtein(this.pdbFiles{this.gui.pdbFiles.Value(1)});
                end
                
            end
            
        end
        
        function view(this, state, trigger, varargin)
            
            switch state
                case 'empty'
                    switch trigger
                        case 'initialize'
                            plotEdProfile(this.gui.ax2);
                        case 'load-data'
                            displayDataFiles();
                            upperPlot();
                        otherwise
                            sprintf('State: %s and trigger: %s is not found for the view.', state, trigger);
                    end
                case 'inspect'
                    switch trigger
                        case 'load-data'
                            displayDataFiles();
                            upperPlot();
                        case 'load-pdb'
                            displayPdbFiles();
                        case {'choose-data', 'toggle-fresnel', 'toggle-cal', 'toggle-fit'}
                            upperPlot();
                        case 'basic-info'
                            upperPlot();
                            lowerPlot();
                        case 'layer-table'
                            plotEdProfile(this.gui.ax2);
                            upperPlot();
                        case 'fit'
                            this.gui.showFit.Value = 1;
                            upperPlot();
                            recordFittingResults();
                        case {'quick-fit', 'angle-fit'}
                            this.gui.showFit.Value = 1;
                            upperPlot();
                            recordFittingResults();
                        case 'update-starts'
                            newStarts = varargin{1};
                            this.gui.parametersTable.Data(:, 3) = num2cell(newStarts');
                        case 'update-starts-follow-up'
                            if this.gui.showCal.Value
                                upperPlot();
                                lowerPlot();
                            end
                        case 'parameter-table'
                            theta = varargin{1};
                            phi = varargin{2};
                            upperPlot();
                            lowerPlot(theta, phi);
                            recordProteinDimension(theta, phi);
                        case 'show-protein'
                            if ~isempty(this.protein)
                                this.gui.showEd.Value = ~ this.gui.showProtein.Value;
                                lowerPlot();
                            else
                                this.gui.showProtein.Value = false;
                            end
                        case 'show-ed'
                            this.gui.showProtein.Value = ~ this.gui.showEd.Value;
                            lowerPlot();
                        case 'protein-table'
                            tdata = varargin{1};
                            theta = tdata.theta;
                            phi = tdata.phi;
                            lowerPlot(theta, phi);
                        case 'toggle-protein'
                            proteinOn = varargin{1};
                            toggleProteinPara(proteinOn);
                        case 'toggle-protein-update'
                            upperPlot();
                            lowerPlot();
                    end
                otherwise
                    disp('Case: %s is not found for the view', state);
                
            end
            
            % utility
            
            function displayDataFiles()
                
                files = cell(1, length(this.data));
                for i = 1 : length(this.data)
                    files{i} = this.data{i}.file;
                end
                
                this.gui.dataFiles.String = files;
                
            end
            
            function displayPdbFiles()
                
                n = length(this.pdbFiles);
                files = cell(1, n);
                for i = 1 : n
                    [~, name, extension] = fileparts(this.pdbFiles{i});
                    files{i} = [name, extension];
                end
                
                this.gui.pdbFiles.String = files;
                
            end
            
            function toggleProteinPara(pro_on)
                
                table = this.gui.parametersTable;
                if pro_on
                    new_row_names = {'Density'; 'Insertion'; 'Theta'; 'Phi'; 'Inhomogeneity'};
                    new_data = {0, 5, 3, false, false; 0, 20, 0, false, false; 0, 180, 0, false, false; 0, 360, 0, false, false; 0, 0, 0, true, false};
                    table.RowName = [table.RowName; new_row_names];
                    table.Data = [table.Data; new_data];
                else
                    table.RowName = table.RowName(1 : end - 5);
                    table.Data = table.Data(1 : end - 5, :);
                end
                
            end
            
            % plotting
            
            function upperPlot()
                
                ax = this.gui.ax1;
                
                sel = cell2mat(this.gui.parametersTable.Data(:, end));
                
                if sum(sel)
                    plotLikelihoodOrChi2(ax);
                else
                    plotSelectedData(ax);
                end
                
            end
            
            function lowerPlot(varargin)
                
                ax = this.gui.ax2;
                switch this.gui.showProtein.Value
                    case 0
                        plotEdProfile(ax);
                    case 1
                        if nargin == 2
                            theta_local = varargin{1};
                            phi_local = varargin{2};
                        else
                            theta_local = this.gui.parametersTable.Data{end-1, 3};
                            phi_local = this.gui.parametersTable.Data{end, 3};
                        end
                        plotProtein(ax, theta_local, phi_local);
                end
                
            end
            
            function plotLikelihoodOrChi2(ax)
                
                sel = cell2mat(this.gui.parametersTable.Data(:, end));
                
                switch sum(sel)
                    case 1
                        loc = find(sel);
                        ind = sum(this.layers.fits.all.fitted(1 : loc));
                        dat = this.layers.fits.one{ind};
                        plot(ax, dat.para_range, dat.likelihood, 'o', 'markersize', 8, 'linewidth', 2);
                        hold(ax, 'on');
                        plot(ax, dat.gauss.x, dat.gauss.y, '-', 'linewidth', 2);
                        hold(ax, 'off');
                        legend(ax, {'Likelihood', 'Gassian Fit'});
                        xlabel(ax, dat.para_name, 'fontsize', 14, 'interpreter', 'latex');
                        ylabel(ax, 'Likelihood', 'fontsize', 14, 'interpreter', 'latex');
                    case 2
                        loc = find(sel);
                        ind1 = sum(this.layers.fits.all.fitted(1 : loc(1)));
                        ind2 = sum(this.layers.fits.all.fitted(1 : loc(2)));
                        dat = this.layers.fits.two{ind1, ind2};
                        contourf(ax, dat.para_range_1, dat.para_range_2, dat.likelihood);
                        colorbar(ax);
                        hold(ax, 'on');
                        plot(ax, dat.confidence.contour(1, :), dat.confidence.contour(2, :), 'r-', 'linewidth', 2);
                        hold(ax, 'off');
                        xlabel(ax, dat.para_names{1}, 'fontsize', 14, 'interpreter', 'latex');
                        ylabel(ax, dat.para_names{2}, 'fontsize', 14, 'interpreter', 'latex');
                        legend(ax, 'Joint Likilihood', '95% confidence window');
                end
                
            end
            
            function plotSelectedData(ax)
                
                selected = this.gui.dataFiles.Value;
                legends = cell(1, length(selected));
                
                switch this.gui.normalized.Value
                    case 0
                        j = 1;
                        for i = selected
                            q = this.data{i}.q();
                            ydata = this.data{i}.ref;
                            err = this.data{i}.err;
                            errorbar(ax, q, ydata, err, 'o', 'linewidth', 1.5, 'color', this.colors(j));
                            legends{j} = this.data{i}.file;
                            j = j + 1;
                            hold(ax, 'on');
                        end
                        ylabel(ax, 'Reflectivity', 'fontsize', 14, 'interpreter', 'latex');
                    case 1
                        j = 1;
                        for i = selected
                            q = this.data{i}.q();
                            d = this.data{i}.getFND();
                            ydata = d.ref;
                            err = d.err;
                            errorbar(ax, q, ydata, err, 'o', 'linewidth', 1.5, 'color', this.colors(j));
                            legends{j} = this.data{i}.file;
                            j = j + 1;
                            hold(ax, 'on');
                        end
                        ylabel(ax, 'Fresnel Normalized Reflectivity', 'fontsize', 14, 'interpreter', 'latex');
                end
                
                if this.gui.showCal.Value
                    dat = this.layers.getFNR();
                    switch this.gui.normalized.Value
                        case 0
                            plot(ax, dat.q, dat.ref, '-', 'linewidth', 2);
                        case 1
                            plot(ax, dat.q, dat.fnr, '-', 'linewidth', 2);
                    end
                    legends = [legends, 'Calculated'];
                end
                
                if this.gui.showFit.Value && isfield(this.layers.fits, 'all')
                    dat = this.layers.fits.all;
                    switch this.gui.normalized.Value
                        case 0
                            plot(ax, dat.q, dat.ref_fit, '-', 'linewidth', 2);
                        case 1
                            plot(ax, dat.q, dat.ref_fit_fnr, '-', 'linewidth', 2);
                    end
                    legends = [legends, 'Fit'];
                end
                
                xlabel(ax, '$$ Q_z(\AA^{-1}) $$', 'interpreter', 'latex', 'fontsize', 14);
                legend(ax, legends, 'interpreter', 'none');
                hold(ax, 'off');
                
            end
            
            function plotEdProfile(ax)
                
%                 plot(ax, this.layers.profile.z, this.layers.profile.layerEd, '-k', 'linewidth', 2);
%                 hold(ax, 'on');
%                 plot(ax, this.layers.profile.z, this.layers.profile.ed, '-b', 'linewidth', 2);
%                 legend(ax, {'Layer Structure', 'Smoothed With Roughness'});
%                 xlabel(ax, '$$ z (\AA) $$', 'interpreter', 'latex', 'fontsize', 14);
%                 ylabel(ax, '$$ Electron Density (\AA^{-3}) $$', 'interpreter', 'latex', 'fontsize', 14);
%                 hold(ax, 'off');
                
                this.layers.plotEdProfile(ax);
                
            end
            
            function plotProtein(ax, theta, phi)
                
                sel_emphasize = [692, 1180, 2175];
                if ~isempty(this.protein)
                    this.protein.visualize(ax, theta, phi, sel_emphasize);
                end
                
            end
            
            % recording
            
            function recordFittingResults()
                
%                 text = {};
%                 if isfield(this.layers.fits, 'all')
%                     dat = this.layers.fits.all;
%                     n = length(dat.para_fitted);
%                     text_1 = cell(n+4, 1);
%                     text_1{1} = this.textDivider();
%                     text_1{2} = 'Global fit result: ';
%                     text_1{3} = '';
%                     text_1{4} = ['Chi^2: ', num2str(dat.chi2)];
%                     for i = 1 : n
%                         text_1{i+4} = [dat.para_names_fitted{i}, ': ', num2str(dat.para_fitted(i))];
%                     end
%                     text = [text; text_1];
%                     if isfield(this.layers.fits, 'one') && ~ isempty(this.layers.fits.one)
%                         dat = this.layers.fits.one;
%                         text_2 = cell(n+3, 1);
%                         text_2{1} = '';
%                         text_2{2} = 'Single parameter brute force fit:';
%                         text_2{3} = '';
%                         for i = 1 : n
%                             text_2{i+3} = [dat{i}.para_name, ': ', num2str(dat{i}.gauss.para(2)), ', standard deviation: ', num2str(dat{i}.gauss.para(3))];
%                         end
%                         text = [text; text_2];
%                     end
%                 end
%                 
%                 this.gui.output.String = [text; this.gui.output.String];
                
            end
            
            function recordProteinDimension(theta, phi)
                
                d = this.layers.protein.dimension(theta, phi);
                
                text = cell(6, 1);
                text{1} = this.textDivider();
                text{2} = ['Protein dimension at theta = ', num2str(theta), ', phi = ', num2str(phi), ':'];
                text{3} = ['x: ', num2str(d.x)];
                text{4} = ['y: ', num2str(d.y)];
                text{5} = ['z: ', num2str(d.z)];
                text{6} = '';
                
                this.gui.output.String = [text; this.gui.output.String];
                
            end
            
        end
        
    end
    
    methods(Static)
        
        function c = colors(n)
            
            cs = 'kbrmcgy';
            index = mod(n, length(cs));
            if index == 0
                index = length(cs);
            end
            
            c = cs(index);
            
        end
        
        function s = textDivider()
            
            s = repmat('-', 1, 76);
            
        end
        
    end
    
end