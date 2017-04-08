classdef XeRef < handle
    
    properties
        
        data
        gui
        handles
        
    end
    
    methods
        
        function this = XeRef()
            
            createView();
            createController();
            
            function createView()
                
                initializeView();
                createFilePanel();
                createAxes();
                createRightPanel();
                
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
                    
                    this.gui.showError = uicontrol(this.handles, 'Style','checkbox','String','Show Error','Units','normalized','Visible','on',...
                        'Position',[0.6 0.965 0.1 0.018]);
                    
                    this.gui.likelihoodChi2 = uicontrol(this.handles,'Style','popupmenu','String',{'Likelihood','Chi^2'},'Visible','off',...
                        'Units','normalized',...
                        'Position',[0.572 0.97 0.1 0.018]);
                    
                    this.gui.showCal = uicontrol(this.handles,'Style','checkbox','String','Show Calc.','Units','normalized',...
                        'Position',[0.6 0.437 0.08 0.018]);
                    
                    this.gui.showFit = uicontrol(this.handles,'Style','checkbox','String','Show Fit','Units','normalized',...
                        'Position',[0.54 0.437 0.06 0.018]);
                    
                    ax1 = axes('Parent',this.handles,'Units','normalized','Position',[0.215 0.52 0.45 0.44]);
                    ax1.XLim = [0 10];
                    ax1.YLim = [0 10];
                    ax1.XTick = [0 2 4 6 8 10];
                    ax1.YTick = [0 2 4 6 8 10];
                    ax1.XLabel.String = 'x1';
                    ax1.YLabel.String = 'y1';
                    
                    this.gui.ax1 = ax1;
                    
                    % plot region 2
                    ax2 = axes('Parent',this.handles,'Units','normalized','Position',[0.215 0.08 0.45 0.35]);
                    ax2.XLim = [0 10];
                    ax2.YLim = [0 10];
                    ax2.XTick = [0 2 4 6 8 10];
                    ax2.YTick = [0 2 4 6 8 10];
                    ax2.XLabel.String = 'x2';
                    ax2.YLabel.String = 'y2';
                    
                    this.gui.ax2 = ax2;
                    
                end
                
            end
            
            function createRightPanel()
                
                this.gui.rightPanel = uipanel(this.handles, 'Units','normalized','Position',[0.68 0.02 0.31 0.97]);
                
            end
            
            function createController()
                
                this.gui.loadData.Callback = @(varargin) this.control('load-data');
                this.gui.dataFiles.Callback = @(varargin) this.control('choose-data');
                
            end
            
        end
        
        function control(this, trigger, varargin)
            
            state = guiState();
            
            switch state
                case 'empty'
                    switch trigger
                        case 'load-data'
                            this.model(state, trigger);
                            this.view(state, trigger);
                        otherwise
                            disp('State: %s, trigger: %s is not found for the controller', state, trigger);
                    end
                case 'inspect'
                    switch trigger
                        case 'load-data'
                            this.model(state, trigger);
                            this.view(state, trigger);
                        case 'choose-data'
                            this.view(state, trigger);
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
                            loadNewData();
                    end
                case 'inspect'
                    switch trigger
                        case 'load-data'
                            loadNewData();
                    end
                otherwise
                    sprintf('Case: %s is not found for the view', state);
            end
            
            function loadNewData()
                
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
        
        function view(this, state, trigger, varargin)
            
            switch state
                case 'empty'
                    switch trigger
                        case 'load-data'
                            displayDataFiles();
                            plotSelectedData();
                        otherwise
                            sprintf('State: %s and trigger: %s is not found for the view.', state, trigger);
                    end
                case 'inspect'
                    switch trigger
                        case 'load-data'
                            displayDataFiles();
                            plotSelectedData();
                        case 'choose-data'
                            plotSelectedData();
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
            
            function plotSelectedData()
                
                selected = this.gui.dataFiles.Value;
                
                for i = selected
                    errorbar(this.gui.ax1, this.data{i}.rawdata.q, this.data{i}.rawdata.ref, this.data{i}.rawdata.err, 'ob', 'linewidth', 2);
                    hold(this.gui.ax1, 'on');
                end
                hold(this.gui.ax1, 'off');
                xlabel(this.gui.ax1, 'Qz');
                ylabel(this.gui.ax1, 'Fresnel Normalized');
                
            end
            
        end
        
    end
    
end