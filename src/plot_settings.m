
function plot_settings()
    % Set general plot settings for the current axes (gca)
    ax = gca;  % Get the current axes handle
    set(ax, 'FontSize', 30, ...                 % Set font size
            'Box', 'on', ...                    % Turn on the box around the plot
            'TickLabelInterpreter', 'latex');   % Set interpreter for tick labels to LaTeX
    
    % Set text interpreters for axis labels and title
    ax.LineWidth = 1;
    ax.XLabel.Interpreter = 'latex';  % Set X-axis label to use LaTeX
    ax.YLabel.Interpreter = 'latex';  % Set Y-axis label to use LaTeX
    if isprop(ax, 'ZLabel')           % Check if Z-axis label exists (e.g., in 3D plots)
        ax.ZLabel.Interpreter = 'latex';  % Set Z-axis label to use LaTeX
    end
    ax.Title.Interpreter = 'latex';   % Set title to use LaTeX

    legendObj = findobj(gcf, 'Type', 'Legend'); % Find the legend object in the current figure
    if ~isempty(legendObj)
        set(legendObj, 'Interpreter', 'latex'); % Set legend interpreter to LaTeX
        set(legendObj, 'FontSize', 26);
    end

    lineObjs = findobj(ax, 'Type', 'Line'); % Find all line objects
    set(lineObjs, 'LineWidth', 1.5);          % Set line width to 2
end
