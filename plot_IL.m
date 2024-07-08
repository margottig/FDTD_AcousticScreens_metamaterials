function plot_IL(metalens, h3_matrix, href)

    % Create the plot
    dt = .01 / 341 / 2;
    h = href;
    n = length(h);
    f = (0:n-1)/(n*dt);  % Frequency vector

    for ii = 1:size(h3_matrix, 1)
        figure(7);
        semilogx(f, 20*log10(abs(fft(href)) ./ abs(fft(h3_matrix(ii, :)))), 'LineWidth', 1);
        hold on;
    end

    set(gca, 'FontSize', 10);
    title("Insertion Loss of the 3 materialens")
    xlabel('f(Hz)', 'FontSize', 15);
    ylabel('IL(dB)', 'FontSize', 15);
    axis([0 5000 0 40]);

    % Add the legend
    legend(metalens, 'Interpreter', 'none');
    

    fcent = [125 250 500 1000 2000 4000];
    for ii = 1:length(fcent)
        Iref(ii) = sum(abs(fft(href)).^2 .* (f < fcent(ii) * sqrt(2)) .* (f > fcent(ii) / sqrt(2)));
        I(1, ii) = sum(abs(fft(h3_matrix(1, :))).^2 .* (f < fcent(ii) * sqrt(2)) .* (f > fcent(ii) / sqrt(2)));
    end

    figure(7);
    hold on;
    semilogx(fcent, 10*log10(Iref ./ I), 'ok'); % ratio of reference intensity (Iref) to measured intensity (I)

end
