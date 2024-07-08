function plot_received_signals(time, h, h2, h3)
% Normalize the signals
maxx = max([abs(h), abs(h2), abs(h3)]);
h = h ./ maxx;
h2 = h2 ./ maxx;
h3 = h3 ./ maxx;

% Create figure
fig = figure('Color', [1, 1, 1]);
hold on;

% Plot each signal with different colors and shapes
plot(time, h, 'r-', 'DisplayName', 'Receiver 1'); % Red line for Receiver 1
plot(time, h2, 'b--', 'DisplayName', 'Receiver 2'); % Blue dashed line for Receiver 2
plot(time, h3, 'g:', 'DisplayName', 'Receiver 3', 'LineWidth', 2); % Green dotted line for Receiver 3, thicker

% Add titles, labels, and legend
title('Received Signal'); 
xlabel('Time (s)');
ylabel('Normalized Signal');
legend;
grid on;

% Set axes properties
set(gca, 'YLim', [-1.2, 1.2]);
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);

hold off;
end
