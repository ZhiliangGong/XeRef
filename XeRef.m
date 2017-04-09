classdef XeRef < handle
    
    properties
        
        data
        layers
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
                createLayersTable();
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
                    
                    this.handles = figure('Visible','off','Name','XeRay','NumberTitle','off','Units','pixels', 'Position', screenSize, 'Resize', 'on');
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
                    
                    this.gui.loadPdb = uicontrol(panel,'Style','pushbutton','String','Load','Units','normalized',...
                        'Position',[0.035 0.015 0.3 0.032]);
                    
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
                    
                end
            
                function createRightPanel()
                    
                    this.gui.rightPanel = uipanel(this.handles, 'Units','normalized','Position',[0.68 0.02 0.31 0.97]);
                    
                end
                
                function createBasicInfoTalbe()
                    
                    rightPanel = this.gui.rightPanel;
                    
                    rowName = {'Beam Energy (keV)'; 'Cut Qz Below'};
                    colName = {};
                    columnFormat = {'numeric'};
                    columnWidth = {120};
                    tableData = {10; 0.026};
                    
                    this.gui.basicInfoTable = uitable(rightPanel, 'Data', tableData, 'ColumnName', colName, ...
                        'ColumnFormat', columnFormat, 'ColumnEditable', true, 'Units','normalized', ...
                        'ColumnWidth',columnWidth,'RowName',rowName, 'RowStriping','off',...
                        'Position', [0.025 0.925 0.935 0.06], 'TooltipString', 'Press enter to update value.');
                    
                end
                
                function createLayersTable()
                    
                    rightPanel = this.gui.rightPanel;
                    
                    rowName = {'Top', 'Layer-2', 'Layer-1', 'Bottom'};
                    colName = {'Electron Density', 'Thickness (A)', 'Delete'};
                    colFormat = {'numeric', 'numeric', 'logical'};
                    colWidth = {140, 90, 40};
                    tableData = { 0, Inf, false; 0.26, 12, false; 0.44, 10, false; 0.334, Inf, false};
                    
                    this.gui.layerText = uicontrol(rightPanel,'Style','text','String','Layer Structure:','Units','normalized','HorizontalAlignment','left',...
                        'Position',[0.025 0.895 0.8 0.025]);
                    
                    this.gui.layerTable = uitable(rightPanel,'Data', tableData,'ColumnName', colName,...
                        'ColumnFormat', colFormat,'ColumnEditable', true(1, 6), 'Units', 'normalized',...
                        'ColumnWidth',colWidth,'RowName',rowName,'RowStriping','off',...
                        'Position', [0.025 0.7 0.935 0.195]);
                    
                    this.gui.addLayer = uicontrol(rightPanel,'Style','pushbutton','String', 'Add', 'Units','normalized', 'Position', [0.725 0.67 0.11 0.03]);
                    
                    this.gui.deleteLayer = uicontrol(rightPanel,'Style','pushbutton','String', 'Delete','Units','normalized', 'Position', [0.84 0.67 0.12 0.03]);
                    
                end
                
                function createParameterTable()
                    
                    rightPanel = this.gui.rightPanel;
                    
                    rowName = {'Qz-Offset', 'Top-ED', 'Layer-2-ED', 'Layer-1-ED', 'Bottom-Ed', 'Layer-2-Thkns', 'Layer-1-Thkns'};
                    colName = {'Min','Max','Start','Fix','Plot'};
                    colFormat = {'numeric','numeric','numeric','logical','logical'};
                    colWidth = {55 55 55 30 30};
                    tableData = {-0.0001, 0.0001, 0, false, false; 0, 0, 0, true, false; 0.2, 0.3, 0.26, false, false;...
                        0.4, 0.45, 0.44, false, false; 0.335, 0.335, 0.335, true, false; 8, 16, 12, false, false;...
                        6, 14, 10, false, false};
                    
                    this.gui.parametersTableTitle = uicontrol(rightPanel,'Style','text','String','Fitting Parameters:',...
                        'Units','normalized','HorizontalAlignment','left', 'Position', [0.025 0.67 0.3 0.025]);
                    
                    this.gui.parametersTable = uitable(rightPanel,'Data', tableData,'ColumnName', colName,...
                        'ColumnFormat', colFormat,'ColumnEditable', [true true true true true],'Units','normalized',...
                        'ColumnWidth',colWidth,'RowName',rowName,'RowStriping','off', 'Position', [0.025 0.37 0.935 0.3]);
                    
                end
                
                function createFittingControls()
                    
                    rightPanel = this.gui.rightPanel;
                    
                    h = 0.34;
                    
                    this.gui.layerTableTitle = uicontrol(rightPanel,'Style','text','String', 'Fitting Control:','Units','normalized','HorizontalAlignment','left',...
                        'Position',[0.025 h 0.8 0.025]);
                    
                    this.gui.loadPara = uicontrol(rightPanel,'Style','pushbutton','String','Load Para','Units','normalized',...
                        'Position',[0.024 h-0.03 0.17 0.03]);
                    
                    this.gui.savePara = uicontrol(rightPanel,'Style','pushbutton','String','Save Para','Units','normalized',...
                        'Position',[0.19 h-0.03 0.17 0.03]);
                    
                    this.gui.stepInput = uicontrol(rightPanel,'Style','edit','String',20,'Units','normalized',...
                        'HorizontalAlignment','left','Position',[0.62 h-0.03 0.1 0.03]);
                    
                    this.gui.stepText = uicontrol(rightPanel,'Style','text','String','Steps','Units','normalized',...
                        'HorizontalAlignment','left','Position', [0.735 h-0.035 0.08 0.03]);
                    
                    this.gui.fitButton = uicontrol(rightPanel,'Style','pushbutton','String','Fit','Units','normalized',...
                        'Position',[0.82 h-0.03 0.15 0.03]);
                    
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
                        'Position',[0.03 0.07 0.935 0.2]);
                    
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
                
                this.gui.loadData.Callback = @(varargin) this.control('load-data');
                this.gui.dataFiles.Callback = @(varargin) this.control('choose-data');
                
                % figure controls
                this.gui.normalized.Callback = @(varargin) this.control('toggle-fresnel');
                this.gui.showCal.Callback = @(varargin) this.control('toggle-cal');
                this.gui.showFit.Callback = @(varargin) this.control('toggle-fit');
                
                % table controls
                this.gui.addLayer.Callback = @(varargin) this.control('add-layer');
                this.gui.deleteLayer.Callback = @(varargin) this.control('delete-layers');
                this.gui.layerTable.CellEditCallback = @(source, eventdata, varargin) this.control('layer-table', eventdata);
                this.gui.parametersTable.CellEditCallback = @(source, eventdata, varargin) this.control('parameter-table', eventdata);
                
                % fitting
                this.gui.fitButton.Callback = @(varargin) this.control('fit');
                this.gui.updateStartButton.Callback = @(varargin) this.control('update-starts');
                
            end
            
        end
        
        function control(this, trigger, varargin)
            
            state = guiState();
            
            switch state
                case 'empty'
                    switch trigger
                        case 'initialize'
                            tableData = readLayerTable();
                            this.model(state, trigger, tableData);
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
                        case {'choose-data', 'toggle-fresnel', 'toggle-cal'}
                            this.view(state, trigger);
                        case 'add-layer'
                            this.view(state, trigger);
                            tableData = readLayerTable();
                            this.model(state, trigger, tableData);
                            this.view(state, 'add-layer-update');
                        case 'layer-table'
                            eventdata = varargin{1};
                            if layerTableEditValid(eventdata)
                                tableData = readLayerTable();
                                this.model(state, trigger, tableData);
                                this.view(state, trigger);
                            end
                        case 'parameter-table'
                            eventdata = varargin{1};
                            if parameterTableEditValid(eventdata)
                                
                            end
                        case 'fit'
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
                            pData = readParameterTable();
                            this.model(state, trigger, pData);
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
            
            function valid = layerTableEditValid(eventdata)
                
                valid = true;
                ind1 = eventdata.Indices(1);
                ind2 = eventdata.Indices(2);
                dat = this.gui.layerTable.Data;
                [m, ~] = size(dat);
                
                switch ind2
                    case 1
                    case 2
                        if ind1 == 1 || ind1 == m
                            this.gui.layerTable.Data{ind1, ind2} = eventdata.PreviousData;
                            valid = false;
                        end
                    case 3
                end
                
            end
            
            function valid = parameterTableEditValid(eventdata)
                
                valid = true;
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
                        elseif newData > table.Data{ind1, 3}
                            table.Data{ind1, 3} = newData;
                        else
                            table.Data{ind1, 4} = false;
                        end
                    case 2
                        if newData < table.Data{ind1, 1}
                            table.Data{ind1, 1} = newData;
                            table.Data{ind1, 3} = newData;
                            table.Data{ind1, 4} = true;
                        elseif newData < table.Data{ind1, 3}
                            table.Data{ind1, 3} = newData;
                        else
                            table.Data{ind1, 4} = false;
                        end
                    case 3
                        if newData < table.Data{ind1, 1}
                            table.Data{ind1, 1} = newData;
                        elseif newData > table.Data{ind1, 2}
                            table.Data{ind1, 2} = newData;
                        end
                    case 4
                    case 5
                end
                
            end
            
            function t = readLayerTable()
                
                dat = this.gui.layerTable.Data;
                t.ed = fliplr(cell2mat(dat(:, 1))');
                t.thickness = fliplr(cell2mat(dat(:, 2))');
                
            end
            
            function t = readParameterTable()
                
                dat = this.gui.parametersTable.Data;
                [m, ~] = size(dat);
                n_layer = (m - 3) / 2;
                
                mat = cell2mat(dat(:, 1:3));
                sel = [1, m : -1 : m - n_layer + 1, m - n_layer : -1 : 2];
                t.p0 = mat(sel, 3)';
                t.lb = mat(sel, 1)';
                t.ub = mat(sel, 2)';
                
            end
            
        end
        
        function model(this, state, trigger, varargin)
            
            switch state
                case 'empty'
                    switch trigger
                        case 'initialize'
                            dat = varargin{1};
                            this.layers = RefLayers(10, dat.ed, dat.thickness);
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
                        case 'layer-table'
                            dat = varargin{1};
                            this.layers = RefLayers(10, dat.ed, dat.thickness);
                        case 'fit'
                            refData = varargin{1};
                            pData = varargin{2};
                            steps = varargin{3};
                            this.layers.fitData(refData, pData.p0, pData.lb, pData.ub, steps);
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
            
        end
        
        function view(this, state, trigger, varargin)
            
            switch state
                case 'empty'
                    switch trigger
                        case 'initialize'
                            plotEdProfile(this.gui.ax2);
                        case 'load-data'
                            displayDataFiles();
                            plotSelectedData(this.gui.ax1);
                        otherwise
                            sprintf('State: %s and trigger: %s is not found for the view.', state, trigger);
                    end
                case 'inspect'
                    switch trigger
                        case 'load-data'
                            displayDataFiles();
                            plotSelectedData(this.gui.ax1);
                        case {'choose-data', 'toggle-fresnel', 'toggle-cal', 'toggle-fit'}
                            plotSelectedData(this.gui.ax1);
                        case 'layer-table'
                            plotEdProfile(this.gui.ax2);
                            plotSelectedData(this.gui.ax1);
                        case 'fit'
                            this.gui.showFit.Value = 1;
                            plotSelectedData(this.gui.ax1);
                        case 'update-starts'
                            newStarts = varargin{1};
                            updateStarts(newStarts);
                        case 'update-starts-follow-up'
                            if this.gui.showCal.Value
                                plotSelectedData(this.gui.ax1);
                            end
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
            
            function updateStarts(starts)
                
                n_layer = (length(starts) + 1) / 2;
                starts = starts([1, end : -1 : end - n_layer + 1, end - n_layer : -1 : 2]);
                this.gui.parametersTable.Data(:, 3) = num2cell(starts');
                
            end
            
            % plotting
            
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
                
                if this.gui.showFit.Value
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
                
                plot(ax, this.layers.profile.z, this.layers.profile.layerEd, '-k', 'linewidth', 2);
                hold(ax, 'on');
                plot(ax, this.layers.profile.z, this.layers.profile.ed, '-b', 'linewidth', 2);
                legend(ax, {'Layer Structure', 'Smoothed With Roughness'});
                xlabel('$$ z (\AA) $$', 'interpreter', 'latex', 'fontsize', 14);
                ylabel('$$ Electron Density (\AA^{-3}) $$', 'interpreter', 'latex', 'fontsize', 14);
                hold(ax, 'off');
                
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
        
    end
    
end