function [th_max,h] = perfprof(A,th_max,colors,strings,linewidth)

%PERFPROF Produce performance profile.
% [th_max, h] = PERFPROF(A,th_max) produces a
% peformance profile for the data in the M-by-N matrix A,
% where A(i,j) > 0 measures the performance of the j’th solver
% on the i-th problem, with smaller values of A(i,j) denoting
% "better". For each solver theta is plotted against the
% probability that the solver is within a factor theta of
% the best solver over all problems, for theta on the interval
% [1, th_max].
% Set A(i,j) = NaN if solver j failed to solve problem i.
% TH_MAX defaults to the smallest value of theta for which
% all probabilities are 1 (modulo any NaN entries of A).
% h is a vector of handles to the lines with h(j)
% corresponding to the j-th solver.

    [m,n] = size(A);
    minA = min(A, [], 2);
    if nargin < 2
        th_max = max(max(A, [], 2) ./ minA);
    end

    % Use definition at some points theta
    n_intervals = 20;
    p_step = 1; % increase to make plot look smoother.
    theta = linspace(1, th_max, n_intervals);
    for j = 1:n
        T = zeros(1, n_intervals);
        for k = 1:n_intervals
            T(k) = sum(A(:,j) <= theta(k)*minA)/m;
        end
        plot(theta(1:p_step:end), T(1:p_step:end), strings{j},...
             'Color', colors(j,:), 'Linewidth', linewidth);
        hold on
    end
    hold off
    axis([1, th_max, 0, 1.0]);

end