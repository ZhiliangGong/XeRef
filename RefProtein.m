classdef RefProtein < handle
    
    properties
        
        file
        pdb
        
        gridSize = 0.5
        
    end
    
    methods
        
        function this = RefProtein(file)
            
            if nargin == 0
                [filename, pathname] = uigetfile('*', 'Select the pdb file.');
                file = fullfile(pathname, filename);
            end
            
            this.file = file;
            text = fileread(file);
            this.pdb = this.readPdbFile(text);
            
            % get the mass and electrons of each atom as a vector
            n = length(this.pdb.x);
            this.pdb.mass = zeros(1, n);
            this.pdb.electron = zeros(1, n);
            this.pdb.radius = zeros(1, n);
            [massTable, electronTable] = this.getPeriodicTable();
            radiusTable = this.getAtomRadiusTable();
            for i = 1 : n
                this.pdb.mass(i) = massTable(this.pdb.atoms(i));
                this.pdb.electron(i) = electronTable(this.pdb.atoms(i));
                this.pdb.radius(i) = radiusTable(this.pdb.atoms(i));
            end
            
        end
        
        function [ed, thickness, area] = getEdProfile(this, theta, phi)
            
            if nargin == 1
                theta = 0;
                phi = 0;
            end
            
            [ed, thickness, area] = this.calculateEdProfile(this.pdb.x, this.pdb.y, this.pdb.z, this.pdb.radius, this.pdb.electron, theta, phi, this.gridSize);
            
        end
        
        % utility
        
        function total = totalElectron(this)
            
            total = sum(this.pdb.electron);
            
        end
        
        function mw = molecularWeight(this)
            
            mw = sum(this.pdb.mass);
            
        end
        
        function grids = generateGrid(this)
            
            pos_bottom = min(positions - repmat(this.pdb.radius, 3, 1), [], 2) - 1e-10;
            pos_top = max(positions + repmat(this.pdb.radius, 3, 1), [], 2) + 1e-10;
            sizes = pos_top - pos_bottom;
            dimensions = ceil(sizes / this.gridSize);
            
            grids = cell(1, 3);
            
            for i = 1 : 3
                grids{i} = pos_bottom(i, :) + this.gridSize * (0 : dimensions(i));
            end
            
        end
        
        function d = dimension(this, theta, phi)
            
            if nargin == 1
                d.x = range(this.pdb.x);
                d.y = range(this.pdb.y);
                d.z = range(this.pdb.z);
            else
               [x, y, z] = this.rotateThetaPhi(this.pdb.x, this.pdb.y, this.pdb.z, theta, phi);
               d.x = range(x);
               d.y = range(y);
               d.z = range(z);
            end
            
        end
        
        function h = height(this, theta, phi)
            
            if nargin == 1
                theta = 0;
                phi = 0;
            end
            
            d = this.dimension(theta, phi);
            h = d.z;
            
        end
        
        function d = getAreaEdProfile(this, theta, phi)
            
            if nargin == 1
                theta = 0;
                phi = 0;
            end
            
            [ed, thick, area] = this.getEdProfile(theta, phi);
            
            d.z = (0 : length(area) - 1) * thick(1);
            d.ed = ed;
            d.area = area;
            
        end
        
        function a = getProjectionArea(this, theta, phi)
            
            if nargin == 1
                theta = 0;
                phi = 0;
            end
            
            [~, ~, area] = this.getEdProfile(theta, phi);
            
            a = max(area);
            
        end
        
        % plot
        
        function visualize(this, ax, theta, phi, sel_emphasize)
            
            if isempty(ax)
                ax = gca;
            end
            
            if isempty(theta)
                theta = 0;
            end
            
            if isempty(phi)
                phi = 0;
            end
            
            sel = this.pdb.atoms == 'C';
            
            [xdata, ydata, zdata] = this.rotateThetaPhi(this.pdb.x(sel), this.pdb.y(sel), this.pdb.z(sel), theta, phi);
            
            plot3(ax, xdata, ydata, zdata, '.', 'linewidth', 2, 'MarkerFaceColor', [1 1 1] * 0.5, 'color', [1 1 1] * 0.5);
            
            if ~ isempty(sel_emphasize)
                [xdata, ydata, zdata] = this.rotateThetaPhi(this.pdb.x(sel_emphasize), this.pdb.y(sel_emphasize), this.pdb.z(sel_emphasize), theta, phi);
                hold(ax, 'on');
                plot3(ax, xdata, ydata, zdata, 's', 'linewidth', 2, 'MarkerFaceColor', 'm', 'color', 'r', 'markersize', 8);
                hold(ax, 'off');
            end
            
            xlabel('$$ y (\AA) $$', 'fontsize', 14, 'interpreter', 'latex');
            ylabel('$$ x (\AA) $$', 'fontsize', 14, 'interpreter', 'latex');
            zlabel('$$ z (\AA) $$', 'fontsize', 14, 'interpreter', 'latex');
            
        end
        
        function contourHeight(this)
            
            this.mapHeight([], [], [], 'contour');
            
        end
        
        function d3Height(this)
            
            this.mapHeight();
            
        end
        
        function ax = mapHeight(this, theta, phi, ax, type)
            
            if nargin == 1 || isempty(theta) || isempty(phi)
                theta = 0 : 1 : 180;
                phi = 0 : 1 : 359;
            end
            
            if nargin < 4 || isempty(ax)
                figure;
                ax = gca;
            end
            
            if nargin < 5 || isempty(type)
                type = 'surf';
            end
            
            [Phi, Theta] = meshgrid(phi, theta);
            m = length(theta);
            n = length(phi);
            heights = zeros(m, n);
            
            for i = 1 : m
                for j = 1 : n
                    heights(i, j) = this.height(theta(i), phi(j));
                end
            end
            
            switch type
                case 'surf'
                    surf(ax, Phi, Theta, heights);
                    xlabel(ax, '\phi (deg.)', 'fontsize', 24);
                    ylabel(ax, '\theta (deg.)', 'fontsize', 24);
                    zlabel(ax, 'Height (\AA)', 'fontsize', 24, 'interpreter', 'latex');
                    set(ax, 'fontsize', 24, 'ylim', [0 180], 'ytick', (0 : 60 : 180), 'xlim', [0 360], 'xtick', (0 : 60 : 360));
                case 'contour'
                    contourf(ax, Phi, Theta, heights);
                    xlabel(ax, '\phi (deg.)', 'fontsize', 24);
                    ylabel(ax, '\theta (deg.)', 'fontsize', 24);
                    set(ax, 'fontsize', 24, 'ylim', [0 180], 'ytick', (0 : 60 : 180), 'xlim', [0 360], 'xtick', (0 : 60 : 360));
            end
            
        end
        
        function plotEdProfile(this, theta, phi, ax)
            
            if nargin == 1
                theta = 0;
                phi = 0;
            end
            
            if nargin < 4
                figure;
                ax = gca;
            end
            
            [ed, thick, ~] = this.getEdProfile(theta, phi);
            
            z_pos = (0 : length(ed) - 1) * thick(1);
            
            hold(ax, 'on');
            plot(ax, z_pos, ed, 'linewidth', 2);
            xlabel('z (\AA)', 'interpreter', 'latex', 'fontsize', 16);
            ylabel('$$ Electron Density (\AA^{-3}) $$', 'interpreter', 'latex', 'fontsize', 16);
            set(ax, 'fontsize', 14);
            hold(ax, 'off');
            
        end
        
        function plotAreaProfile(this, theta, phi, ax)
            
            if nargin == 1
                theta = 0;
                phi = 0;
            end
            
            if nargin < 4
                figure;
                ax = gca;
            end
            
            [~, thick, area] = this.getEdProfile(theta, phi);
            
            z_pos = (0 : length(area) - 1) * thick(1);
            
            hold(ax, 'on');
            plot(ax, z_pos, area, 'linewidth', 2);
            xlabel('z (\AA)', 'interpreter', 'latex', 'fontsize', 16);
            ylabel('$$ Protein Area (\AA^{2}) $$', 'interpreter', 'latex', 'fontsize', 16);
            set(ax, 'fontsize', 14);
            hold(ax, 'off');
            
        end
        
    end
    
    methods(Static)
        
        function [AtomicMassTable, ElectronTable] = getPeriodicTable()

            symbols = 'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U';
            atomicNumber = (1:92);
            masses = [1.0079,4.0026,6.941,9.0122,10.81,12.011,14.007,15.999,18.998,20.18,22.99,24.305,26.982,28.085,30.974,32.066,35.453,39.948,39.098,40.078,44.956,47.867,50.941,51.996,54.938,55.845,58.933,58.693,63.546,65.39,69.723,72.61,74.922,78.96,79.904,83.8,85.468,87.62,88.906,91.224,92.906,95.94,98,101.07,102.91,106.42,107.87,112.41,114.82,118.71,121.76,127.6,126.9,131.29,132.91,137.33,138.91,140.12,140.91,144.24,145,150.36,151.96,157.25,158.93,162.5,164.93,167.26,168.93,173.04,174.97,178.49,180.95,183.84,186.21,190.23,192.22,195.08,196.97,200.59,204.38,207.2,208.98,209,210,222,223,226,227,232.04,231.04,238.03];

            symbol = regexp(symbols, ' ', 'split');
            ElectronTable = containers.Map;
            AtomicMassTable = containers.Map;

            for i = 1:length(symbol)
                ElectronTable(symbol{i}) = atomicNumber(i);
                AtomicMassTable(symbol{i}) = masses(i);
            end

        end
        
        function radius = getAtomRadiusTable()
            
            radius = containers.Map;
            radius('C') = 1.7;
            radius('N') = 1.55;
            radius('H') = 1.1;
            radius('O') = 1.52;
            radius('S') = 1.8;
            radius('P') = 1.95;
            
        end
        
        function [x, y, z] = rotateThetaPhi(x, y, z, theta, phi)
            
            theta = theta * pi / 180;
            phi = phi * pi / 180;
            
            positions = [1, 0, 0; 0, cos(theta), -sin(theta); 0, sin(theta), cos(theta)]...
                * [cos(phi), -sin(phi), 0; sin(phi), cos(phi), 0; 0, 0, 1]...
                * [x; y; z];
            
            x = positions(1, :);
            y = positions(2, :);
            z = positions(3, :);
            
        end
        
        function [ed, thickness, area] = calculateEdProfile(x, y, z, radius, electron, theta, phi, gs)
            
            N = length(x);
            
            theta = theta * pi / 180;
            phi = phi * pi / 180;
            
            positions = [1, 0, 0; 0, cos(theta), -sin(theta); 0, sin(theta), cos(theta)]...
                * [cos(phi), -sin(phi), 0; sin(phi), cos(phi), 0; 0, 0, 1]...
                * [x; y; z];
            
            x = positions(1, :);
            y = positions(2, :);
            z = positions(3, :);
            
            xtop = max(x + radius) + 1e-10;
            xbot = min(x - radius) - 1e-10;
            ytop = max(y + radius) + 1e-10;
            ybot = min(y - radius) - 1e-10;
            ztop = max(z + radius) + 1e-10;
            zbot = min(z - radius) - 1e-10;
            
            xn = ceil((xtop - xbot) / gs);
            yn = ceil((ytop - ybot) / gs);
            zn = ceil((ztop - zbot) / gs);
            
            x_grid_pos = xbot + gs * ((-1 : xn) + 0.5);
            y_grid_pos = ybot + gs * ((-1 : yn) + 0.5);
            z_grid_pos = zbot + gs * ((-1 : zn) + 0.5);
            
            grid_x = repmat(x_grid_pos, N, 1);
            [x_top_indices(:, 1), x_top_indices(:, 2)] = find(grid_x >= (x + radius)' & grid_x < (x + radius + gs)');
            [x_bot_indices(:, 1), x_bot_indices(:, 2)] = find(grid_x <= (x - radius)' & grid_x > (x - radius - gs)');
            
            grid_y = repmat(y_grid_pos, N, 1);
            [y_top_indices(:, 1), y_top_indices(:, 2)] = find(grid_y >= (y + radius)' & grid_y < (y + radius + gs)');
            [y_bot_indices(:, 1), y_bot_indices(:, 2)] = find(grid_y <= (y - radius)' & grid_y > (y - radius - gs)');
            
            grid_z = repmat(z_grid_pos, N, 1);
            [z_top_indices(:, 1), z_top_indices(:, 2)] = find(grid_z >= (z + radius)' & grid_z < (z + radius + gs)');
            [z_bot_indices(:, 1), z_bot_indices(:, 2)] = find(grid_z <= (z - radius)' & grid_z > (z - radius - gs)');
            
            x_top_indices = sortrows(x_top_indices);
            y_top_indices = sortrows(y_top_indices);
            z_top_indices = sortrows(z_top_indices);
            x_bot_indices = sortrows(x_bot_indices);
            y_bot_indices = sortrows(y_bot_indices);
            z_bot_indices = sortrows(z_bot_indices);
            
            [xf_grid, yf_grid, zf_grid] = meshgrid(x_grid_pos, y_grid_pos, z_grid_pos);
            
            xf_grid = permute(xf_grid, [2, 1, 3]);
            yf_grid = permute(yf_grid, [2, 1, 3]);
            zf_grid = permute(zf_grid, [2, 1, 3]);
            
            Elec_grid = zeros(xn + 2, yn + 2, zn + 2);
            for k = 1 : N
                
                Elec_frac = zeros(x_top_indices(k, 2)-x_bot_indices(k, 2)+1, y_top_indices(k, 2)-y_bot_indices(k, 2)+1, z_top_indices(k, 2)-z_bot_indices(k, 2)+1);
                
                atomdist = (xf_grid(x_bot_indices(k, 2):x_top_indices(k, 2),y_bot_indices(k, 2):y_top_indices(k, 2), z_bot_indices(k, 2):z_top_indices(k, 2)) - x(k)).^2 ...
                    + (yf_grid(x_bot_indices(k, 2):x_top_indices(k, 2),y_bot_indices(k, 2):y_top_indices(k, 2), z_bot_indices(k, 2):z_top_indices(k, 2)) - y(k)).^2 ...
                    + (zf_grid(x_bot_indices(k, 2):x_top_indices(k, 2),y_bot_indices(k, 2):y_top_indices(k, 2), z_bot_indices(k, 2):z_top_indices(k, 2)) - z(k)).^2;
                
                ind_list = (atomdist <= radius(k)^2);
                
                Elec_frac(ind_list) = (electron(k) * gs^3)/(4 / 3 * pi * radius(k)^3);
                
                Elec_grid(x_bot_indices(k, 2):x_top_indices(k, 2),y_bot_indices(k, 2):y_top_indices(k, 2), z_bot_indices(k, 2):z_top_indices(k, 2))...
                    = Elec_grid(x_bot_indices(k, 2):x_top_indices(k, 2),y_bot_indices(k, 2):y_top_indices(k, 2), z_bot_indices(k, 2):z_top_indices(k, 2)) + Elec_frac;
                
            end
            
            thickness = ones(1, zn + 2) * gs;
            area = squeeze(sum(sum(Elec_grid > 0, 1), 2))' * gs^2;
            area(area == 0) = 1; % to avoid division by zero
            ed = squeeze(sum(sum(Elec_grid, 1), 2))' ./ area / gs;
            
            sel_ed = ed > 0;
            ed = ed(sel_ed);
            area = area(sel_ed);
            
        end
        
        function pdb = readPdbFile(text)
            
            raw = textscan(text, '%s %d %s %*s %*s %d %f %f %f %*f %*f %*s', 'HeaderLines', 1);
            
            labels = raw{1};
            anisous = raw{2};
            
            selects = (strcmpi('atom', labels) | strcmpi('hetatm', labels)) & anisous;
            pdb.resid = raw{4}(selects)';
            pdb.x = raw{5}(selects)';
            pdb.y = raw{6}(selects)';
            pdb.z = raw{7}(selects)';
            
            atomNames = raw{3}(selects)';
            atoms = repmat(' ', 1, length(atomNames));
            
            for i = 1 : length(atoms)
                atoms(i) = atomNames{i}(1);
            end
            
            pdb.atoms = atoms;
            
        end
        
    end
    
end