format compact

addpath('external');
addpath('include');
warning('off');

n = 10;

initialize_tests

ids_min = 1;
[~, ids_max] = mygallery();

p = 16;
pp = 2 * p;

produce_results = true;
plot_out = true;
png_out = true;
tikz_out = true;
compute_condest = true;

if produce_results

    n_mats = ids_max - ids_min + 1;

    bwbound_expm = zeros(1, n_mats);
    fwderr_expm = zeros(1, n_mats);
    fwderr_expotf = zeros(1, n_mats);
    fwderr_expp = zeros(1, n_mats);
    fwderr_expfp = zeros(1, n_mats);
    fwderr_expt = zeros(1, n_mats);
    fwderr_expft = zeros(1, n_mats);
    condest1u = zeros(1, n_mats);
    hermitian_select = zeros(1, n_mats);

    group_ids = cell(n_mats, 1);

    rng(1);

    for i = ids_min : ids_max

        fprintf('\n* Matrix id: %d\n', i);

        A = mygallery(i,n);
        if ishermitian(A)
            hermitian_select(i) = 1;
        end

        mp.Digits(p);
        fprintf('Processing MATLAB...\n');
        X_expm = double(expm(A));
        fprintf('Processing expm_otf...\n');
        [X_expotf, m_otf, s_otf, Cp, b] = exptayotf(A);
        fprintf('Processing expp...\n');
        [X_expp, s_expp, mxm_expp] = expm_mp(double(A),...
                                             'algorithm', 'complexschur',...
                                             'abserr', false,...
                                             'approx', 'diagonalcheap',...
                                             'maxscaling', 50,...
                                             'maxdegree', 50);
        fprintf('Processing expfp...\n');
        [X_expfp, s_expfp, mxm_expfp] = expm_mp(double(A),...
                                                'algorithm', 'transfree',...
                                                'abserr', false,...
                                                'approx', 'diagonalcheap',...
                                                'maxscaling', 50,...
                                                'maxdegree', 50);
        fprintf('Processing expt...\n');
        [X_expt, s_expt, mxm_expt] = expm_mp(double(A),...
                                             'algorithm', 'complexschur',...
                                             'abserr', false,...
                                             'approx', 'taylor',...
                                             'maxscaling', 50,...
                                             'maxdegree', 150);
        fprintf('Processing expft...\n');
        [X_expft, s_expft, mxm_expft] = expm_mp(double(A),...
                                                'algorithm', 'transfree',...
                                                'abserr', false,...
                                                'approx', 'taylor',...
                                                'maxscaling', 50,...
                                                'maxdegree', 150);

        mp.Digits(pp);
        exact = expm(mp(A));

        fwderr_expm(i) = double(norm(X_expm - exact, 1) / norm(exact, 1));
        fwderr_expotf(i) = double(norm(X_expotf - exact, 1) / norm(exact, 1));
        fwderr_expp(i) = double(norm(X_expp - exact, 1) / norm(exact, 1));
        fwderr_expfp(i) = double(norm(X_expfp - exact, 1) / norm(exact, 1));
        fwderr_expt(i) = double(norm(X_expt - exact, 1) / norm(exact, 1));
        fwderr_expft(i) = double(norm(X_expft - exact, 1) / norm(exact, 1));

        fprintf('Processing condest...\n');
        if isnan(fwderr_expp(i))
            i
        end
        if (compute_condest)
            condest1u(i) = expm_cond(A) * eps();
        end

    end
    save('test_expm_data_double.mat',...
         'hermitian_select', 'condest1u',...
         'fwderr_expm', 'fwderr_expotf',...
         'fwderr_expt', 'fwderr_expft',...
         'fwderr_expp', 'fwderr_expfp');
else
    load test_expm_data_double.mat
end





%% Plots
print_legend = false;

Tcolors = [color_expm; color_expotf;...
           color_expt; color_expft; color_expp; color_expfp];
Tstyles = {ls_expm; ls_expotf;...
           ls_expt; ls_expft; ls_expp; ls_expfp};
Tmarkers = {marker_expm; marker_expotf;...
            marker_expt; marker_expft; marker_expp; marker_expfp};

legend_perm = [1 2 3 4 5 6];

indices1 = find(~hermitian_select);
[~, perm1] = sort(condest1u(indices1), 'descend');
indices2 = find(hermitian_select);
[~, perm2] = sort(condest1u(indices2), 'descend');

% Plot performance profile
ymin = 1e-18;
ymax = 1e-3;

produce_perfprof(fwderr_expm(indices1), fwderr_expotf(indices1),...
                 fwderr_expt(indices1), fwderr_expft(indices1),...
                 fwderr_expp(indices1), fwderr_expfp(indices1),...
                 16, Tcolors, Tstyles, print_legend, legend_perm,...
                 png_out, tikz_out, 'noherm_double');

produce_error_plot(condest1u(indices1), perm1, ...
                   fwderr_expm(indices1), fwderr_expotf(indices1),...
                   fwderr_expt(indices1), fwderr_expft(indices1),...
                   fwderr_expp(indices1), fwderr_expfp(indices1),...
                   16, ymin, ymax,...
                   Tcolors, color_cond, Tmarkers, ls_cond, msize, lw,...
                   true, print_legend, legend_perm,...
                   'noherm', png_out, tikz_out, 'noherm_double');

produce_error_plot(condest1u(indices2), perm2, ...
                   fwderr_expm(indices2), fwderr_expotf(indices2),...
                   fwderr_expt(indices2), fwderr_expft(indices2),...
                   fwderr_expp(indices2), fwderr_expfp(indices2),...
                   16, ymin, ymax,...
                   Tcolors, color_cond, Tmarkers, ls_cond, msize, lw,...
                   0, false, false,...
                   'herm', png_out, tikz_out, 'herm_double');