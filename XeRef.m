classdef XeRef < handle
    
    properties
        
        data
        gui
        handles
        
    end
    
    methods
        
        function this = XeRef(files)
            
            createView();
            createController();
            
            function createView()
                
                initializeView();
                createFilePanel();
                createAxes();
                createFigureControls();
                createRightPanel();
                createLayersTable();
                createBasicInfoTalbe();
                
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
                        'Position',[0.6 0.915 0.1 0.018]);
                    
                    this.gui.likelihoodChi2 = uicontrol(this.handles,'Style','popupmenu','String',{'Likelihood','Chi^2'},'Visible','off',...
                        'Units','normalized',...
                        'Position',[0.572 0.97 0.1 0.018]);
                    
                    this.gui.showFit = uicontrol(this.handles,'Style','checkbox','String','Show Fit','Units','normalized',...
                        'Position',[0.54 0.437 0.06 0.018]);
                    
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
                    
                    rowName = {'Top', 'Bottom'};
                    colName = {'Electron Density', 'Thickness (A)', 'Delete'};
                    colFormat = {'numeric', 'numeric', 'logical'};
                    colWidth = {140, 90, 40};
                    tableData = { 0, Inf, false; 0.334, Inf, false};
                    
                    this.gui.layerText = uicontrol(rightPanel,'Style','text','String','Layer Structure:','Units','normalized','HorizontalAlignment','left',...
                        'Position',[0.025 0.895 0.8 0.025]);
                    
                    this.gui.layerTable = uitable(rightPanel,'Data', tableData,'ColumnName', colName,...
                        'ColumnFormat', colFormat,'ColumnEditable', true(1, 6), 'Units', 'normalized',...
                        'ColumnWidth',colWidth,'RowName',rowName,'RowStriping','off',...
                        'Position', [0.025 0.745 0.935 0.15]);
                    
                    this.gui.addLayer = uicontrol(rightPanel,'Style','pushbutton','String', 'Add', 'Units','normalized', 'Position', [0.725 0.715 0.11 0.03]);
                    
                    this.gui.deleteLayer = uicontrol(rightPanel,'Style','pushbutton','String', 'Delete','Units','normalized', 'Position', [0.84 0.715 0.12 0.03]);
                    
                end
                
            end
            
            function createController()
                
                this.gui.loadData.Callback = @(varargin) this.control('load-data');
                this.gui.dataFiles.Callback = @(varargin) this.control('choose-data');
                
                % figure controls
                this.gui.normalized.Callback = @(varargin) this.control('toggle-fresnel');
                this.gui.showCal.Callback = @(varargin) this.control('toggle-cal');
                
            end
            
            if nargin == 1
                this.control('load-data', files);
            end
            
        end
        
        function control(this, trigger, varargin)
            
            state = guiState();
            
            switch state
                case 'empty'
                    switch trigger
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
                        otherwise
                            sprintf('State: %s, trigger: %s is not found for the controller', state, trigger)
                    end
                otherwise
                    sprintf('State: %s not found for the controller', state);
            end
            
            function state = guiState()
                
                if isempty(this.gui.dataFiles.String)
                    state = 'empty';
                else
                    state = 'inspect';
                end
                
            end
            
            
            
        end
        
        function model(this, state, trigger, varargin)
            
            switch state
                case 'empty'
                    switch trigger
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
                    end
                otherwise
                    sprintf('Case: %s is not found for the view', state);
            end
            
            function loadNewData(datafiles)
                
                if nargin == 1
                    files = datafiles;
                    newdata = cell(1, length(files));
                    
                    for i = 1 : length(files)
                        newdata{i} = XeRefData(files{i});
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
                            newdata{i} = XeRefData(fullfile(path, files{i}));
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
                        case {'choose-data', 'toggle-fresnel', 'toggle-cal'}
                            plotSelectedData(this.gui.ax1);
                    end
                otherwise
                    disp('Case: %s is not found for the view', state);
                
            end
            
            function displayDataFiles()
                
                files = cell(1, length(this.data));
                for i = 1 : length(this.data)
                    files{i} = this.data{i}.rawdata.file;
                end
                
                this.gui.dataFiles.String = files;
                
            end
            
            function plotSelectedData(ax)
                
                selected = this.gui.dataFiles.Value;
                legends = cell(1, length(selected));
                
                switch this.gui.normalized.Value
                    case 0
                        j = 1;
                        for i = selected
                            q = this.data{i}.layers.getQ();
                            ydata = this.data{i}.layers.data.ref;
                            err = this.data{i}.layers.data.err;
                            errorbar(ax, q, ydata, err, 'o', 'linewidth', 2, 'color', this.colors(j));
                            legends{j} = this.data{i}.rawdata.file;
                            j = j + 1;
                            hold(ax, 'on');
                        end
                        ylabel(ax, 'Reflectivity', 'fontsize', 14);
                    case 1
                        j = 1;
                        for i = selected
                            q = this.data{i}.layers.getQ();
                            d = this.data{i}.layers.getFresnelNormalizedData();
                            ydata = d.ref;
                            err = d.err;
                            errorbar(ax, q, ydata, err, 'o', 'linewidth', 2, 'color', this.colors(j));
                            legends{j} = this.data{i}.rawdata.file;
                            j = j + 1;
                            hold(ax, 'on');
                        end
                        ylabel(ax, 'Fresnel Normalized Reflectivity', 'fontsize', 14);
                end
                
                if this.gui.showCal.Value
                    switch this.gui.normalized.Value
                        case 0
                            j = 1;
                            for i = selected
                                q = this.data{i}.layers.getQ();
                                ydata = this.data{i}.layers.getRef();
                                plot(ax, q, ydata, '-', 'Color', this.colors(j), 'linewidth', 2);
                            end
                        case 1
                            j = 1;
                            for i = selected
                                q = this.data{i}.layers.getQ();
                                ydata = this.data{i}.layers.getFNR();
                                plot(ax, q, ydata, '-', 'Color', this.colors(j), 'linewidth', 2);
                            end
                    end
                end
                
                xlabel(ax, '$$ Q_z(\AA^{-1}) $$', 'interpreter', 'latex', 'fontsize', 14);
                legend(ax, legends, 'interpreter', 'none');
                hold(ax, 'off');
                
            end
            
        end
        
    end
    
    methods(Static)
        
        function color = colors(n)
            
            colors = 'kbrmcgy';
            index = mod(n, length(colors));
            if index == 0
                index = length(colors);
            end
            
            color = colors(index);
            
        end
        
    end
    
end