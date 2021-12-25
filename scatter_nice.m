function scatter_nice(xdata, ydata, varargin)
%% SCATTER_NICE scatter plot wrapper with benefits
% Wrapper for scatter(.) that
% 1) chooses useful defaults for colors and marker size and shape,
% 2) allows for easy group-wise scatter plots, and
% 3) implements randomized-order plotting to prevent the last plotted data to (misleadingly)
%    dominate plot appearance,
% 4) uses nonlinear color scaling to exploit the full color range, and
% 5) implements some automatic simple statistical annotations to be shown (optionally).
%
% The goal is for this to be a one-stop function to creating a useful and visually appealing
% scatter diagram.
%
% Requires the `cbrewer` file exchange function to be on the path, see 
% https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab.
% The file is included here for convenience; copyright remains with the original author.
%
% MANDATORY ARGUMENTS
% xdata: either a vector or (in the case of grouped data) a cell array of vectors.
% ydata: same as xdata.
%
% OPTIONAL POSITIONAL ARGUMENTS [can be omitted, but, if provided, must be provided in that order]
% cdata: color specifications. In the case of just a single group, this can have all the formats 
%   accepted by scatter(.). In the case of group-wise data, this must also be a cell array, with
%   each cell providing the color specification for the corresponding group in a format accepted by
%   scatter(.). If not provided or empty (cdata={}), a qualitative color map with the correct
%   number of groups is used. In the case of a single group and length(cdata)==length(xdata),
%   each datum is colored along a linear color map according to its cdata value. A diverging color 
%   map is chosen in this case. 
%   The two color maps are the 'Set1' and 'RdYlBu' color maps from https://colorbrewer2.org.
% x_label: x axis label.
% y_label: y axis label.
% plot_title: title to put above the plot.
% group_names: cell array with as many cells as xdata and ydata. Each cell must
%    contain the label of the correponding data group, to be used for the plot legend.
% 
% NAME-VALUE ARGUMENTS [optional, must come after positional args, in "..., 'Name', Val, ..." format]
% analyze: if true, add statistical annotation to the diagram. Default: false.
% randomize_order: plot data in randomized order to prevent, e.g. everything being overshadowed by
%    the last group. Default: true.
% balance_colors: if true (default) and linear coloring is used (see above), balance the color map
%    based on the data set such that all colors appear with (roughly) equal frequency.
%
% Eike Petersen, 2020-2021
%

    %% Input handling
    
    if ~iscell(xdata)
        assert(~iscell(ydata));
        xdata = {xdata};
        ydata = {ydata};
    end
    
    p = inputParser;

    % Positional arguments
    addRequired(p, 'xdata', @iscell);
    addRequired(p, 'ydata', @iscell);
    addOptional(p, 'cdata', {});
    addOptional(p, 'x_label', 'x', @ischar);
    addOptional(p, 'y_label', 'y', @ischar);
    addOptional(p, 'plot_title', '', @ischar);
    addOptional(p, 'group_names', {}, @iscell);
    
    % Name-value arguments
    addParameter(p, 'analyze', false, @islogical);
    addParameter(p, 'randomize', true, @islogical);
    addParameter(p, 'balance_colors', true, @islogical);
    
    parse(p, xdata, ydata, varargin{:})
    struct2vars(p.Results);       
    
    if ~isempty(cdata) && length(cdata) == length(xdata{1})
        % Data are colored individually. Pick a diverging color scheme.
        color_data_provided = true;
        cdata = {cdata};
        balcdata = cell(length(xdata), 1);
    else
        % Pick a qualitative color scheme with as many colors as there are groups to plot
        % Will _not_ be used if there is only one group _and_ cdata are provided.
        % (See further below for that case.)
        c_rgb = cbrewer('qual', 'Set1', length(xdata));
        color_data_provided = false;
    end
    
    colors = cell(length(xdata), 1);
    for ii = 1:length(xdata)  % number of groups

        if isempty(cdata)
            % Plot each group in a different color
            colors{ii} = c_rgb(ii, :);
        else
           % Use user-provided colors.
            if color_data_provided
                if balance_colors
                    [balcdata{ii}, cbhticks, cbhticklabels] = balance_cdata(cdata{ii});
                    colors{ii} = balcdata{ii};
                else
                    colors{ii} = cdata{ii};
                end
            else
                colors = cdata{ii};
            end
        end
    end
    
    % choose marker size adaptively, based on the number of points
    N = length(xdata) * length(xdata{1});
    sz = max(45-10*log10(N), 5);
    
    handles = zeros(size(xdata));
    
    if ~randomize
        % Plot data points group-wise. This may lead to the data belonging
        % to the last group hiding all others!
        for ii = 1:length(xdata)  % number of groups
            
            if ~isempty(group_names)
                if sz < 15
                    handles(ii) = scatter(xdata{ii}, ydata{ii}, sz, colors{ii}, 'filled', ...
                        'DisplayName', group_names{ii});
                else
                    handles(ii) = scatter(xdata{ii}, ydata{ii}, sz, colors{ii}, ...
                        'DisplayName', group_names{ii});
                end
            else
                if sz < 15
                    handles(ii) = scatter(xdata{ii}, ydata{ii}, sz, colors{ii}, 'filled');
                else
                    handles(ii) = scatter(xdata{ii}, ydata{ii}, sz, colors{ii});
                end
            end
            hold on;
        end
    else
        % Plot data in a random order to prevent one group - misleadingly -
        % dominating the figure, simply because it was plotted last.
        
        % First, plot one datum from each group to set up legend entries
        for ii = 1:length(xdata)  % number of groups
            
            if color_data_provided
                first_color = colors{ii}(1);
            else
                first_color = colors{ii};
            end
            
            if ~isempty(group_names)
                if sz < 15
                    handles(ii) = scatter(xdata{ii}(1), ydata{ii}(1), sz, first_color, 'filled', ...
                        'DisplayName', group_names{ii});
                else
                    handles(ii) = scatter(xdata{ii}(1), ydata{ii}(1), sz, first_color, ...
                        'DisplayName', group_names{ii});
                end
            else
                if sz < 15
                    handles(ii) = scatter(xdata{ii}(1), ydata{ii}(1), sz, first_color, 'filled');
                else
                    handles(ii) = scatter(xdata{ii}(1), ydata{ii}(1), sz, first_color);
                end
            end
            hold on;
        end
        
        % Now set up vectors containing all data points from all groups
        xarray = [xdata{:}];
        yarray = [ydata{:}];
        if color_data_provided
            c_rgb_array = zeros(size(xarray, 1), 1);
        else
            c_rgb_array = zeros(size(xarray, 1), 3);
        end
        
        % Set up color specification
        kk = 0;
        for ii = 1:length(xdata)
            idces = kk + (1:length(xdata{ii}));
            kk = idces(end);
            if color_data_provided
                % Different colors per datum
                if balance_colors
                    % To make each color appear equally often, use ecdf(color) instead - that's
                    % uniformly distributed.
                    c_rgb_array(idces, :) = balcdata{ii};
                else
                    c_rgb_array(idces, :) = cdata{ii};
                end
            else
                % One color per group
                c_rgb_array(idces, :) = repmat(c_rgb(ii, :), length(xdata{ii}), 1);
            end
        end
        
        % Shuffle vectors
        perm_idces = randperm(length(xarray));
        
        if sz < 15
            scatter(xarray(perm_idces), yarray(perm_idces), sz, c_rgb_array(perm_idces, :), 'filled');
        else
            scatter(xarray(perm_idces), yarray(perm_idces), sz, c_rgb_array(perm_idces, :));
        end
    end
    
    if color_data_provided
        colormap(gca, flipud(cbrewer('div', 'RdYlBu', 11)));
        cbh = colorbar;
        if balance_colors
            % We transformed the color values to ecdf(color) to make all colors appear equally
            % often. Make sure that we display the original data values on the colorbar!
            cbh.Ticks = cbhticks;
            cbh.TickLabels = cbhticklabels;
        end
    end
    
    if ~isempty(x_label)
        xlabel(x_label);
    end
    
    if ~isempty(y_label)
        ylabel(y_label);
    end
    
    if ~isempty(plot_title)
        title(plot_title);
    end
    
    leg = legend(handles, 'interpreter', 'none');
    % Enable clicking on axis to hide lines
    leg.ItemHitFcn = @hileline;
    
    if analyze
        % Statistical annotation
        [m, n] = size(xdata{1});
        if m > n
            all_data_x = cat(1, xdata{:});
            all_data_y = cat(1, ydata{:});
            non_nan_mask = ~any(isnan([all_data_x, all_data_y]'));
        else
            all_data_x = cat(2, xdata{:});
            all_data_y = cat(2, ydata{:});
            non_nan_mask = ~any(isnan([all_data_x; all_data_y]));
        end
        
        r2 = corr(all_data_x(non_nan_mask)', all_data_y(non_nan_mask)', 'Type', 'Pearson');
        rho = corr(all_data_x(non_nan_mask)', all_data_y(non_nan_mask)', 'Type', 'Spearman');

        lsline;

        txt = {['R2: ' num2str(r2)], ...
            ['rho: ', num2str(rho)]};
        text(0.05, 0.95, txt, 'Units', 'normalized', 'VerticalAlignment', 'top');
    end
end

function hileline(src, event)
% This callback toggles the visibility of the line

    if strcmp(event.Peer.Visible,'on')   % If current line is visible
        event.Peer.Visible = 'off';      %   Set the visibility to 'off'
    else                                 % Else
        event.Peer.Visible = 'on';       %   Set the visibility to 'on'
    end
end