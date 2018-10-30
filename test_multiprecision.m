warning('off');

% addpath
addpath('./external')
addpath('./include')
addpath('./matlog')

n = 10;

ids_min = 1;
[~, ids_max] = mygallery();

initialize_tests

% precisions
p_vec = [64, 256, 1024];
ymins = [1e-66, 1e-258, mp('1e-1026')];
ymaxs = [1e-54, 1e-243, mp('1e-1014')];

if (~exist('paused', 'var'))
    paused = false;
end

produce_results = true;
plot_out = true;
png_out = true;
tikz_out = true;

M = ids_max - ids_min + 1;
P = length(p_vec);
condest1 = mp(zeros(1, M));

smax = 100;
mmax_diag = 400;
mmax_tayl = 1000;

if produce_results
    condest1u = mp(zeros(P, M));
    hermitian_select = zeros(1,M);

    fwderr_expmct = mp(zeros(P, M));
    fwderr_expotf = mp(zeros(P, M));

    fwderr_expt = mp(zeros(P, M));
    fwderr_expft = mp(zeros(P, M));
    fwderr_expp = mp(zeros(P, M));
    fwderr_expfp = mp(zeros(P, M));

    main_loop = tic;

    for j = 1:P
        p = p_vec(j);
        pp = floor(p * 2);

        mp.Digits(p);

        s_otf = zeros(1, M);
        s_expt = zeros(1, M);
        s_expft = zeros(1, M);
        s_expp = zeros(1, M);
        s_expfp = zeros(1, M);

        c_expt = zeros(4, M);
        c_expft = zeros(4, M);
        c_expp = zeros(4, M);
        c_expfp = zeros(4, M);

        m_otf = zeros(1, M);
        m_expt = zeros(1, M);
        m_expft = zeros(1, M);
        m_expp = zeros(1, M);
        m_expfp = zeros(1, M);

        n_mats = ids_max - ids_min + 1;

        for k = ids_min : ids_max

            fprintf('\n* Matrix id: %d\n', k);

            mp.Digits(p);
            epsilon = mp('eps');

            A = mygallery(k,n);
            if ishermitian(A)
                hermitian_select(k) = 1;
            end

            % Compute reference solution.
            mp.Digits(pp);
            exact = expm(mp(A));
            mp.Digits(p);

            norm_exact = norm(exact, 1);
            normA = norm(mp(A), 1);

            % Check backward bound of reference solution.
            mp.Digits(pp);
            mp.Digits(p);
            fprintf('Processing expm_mct...\n');
            expmct = expm(mp(A));

            fprintf('Processing expm_otf...\n');
            [expotf, m_otf(k), s_otf(k), Cp, b] = exptayotf(mp(A), eps('mp'));

            fprintf('Processing expt...\n');
            [expt, s_expt(k), m_expt(k), c_expt(:, k)] =...
                expm_mp(mp(A),...
                        'epsilon', eps('mp'),...
                        'algorithm', 'complexschur',...
                        'abserr', false,...
                        'approx', 'taylor',...
                        'maxscaling', smax,...
                        'maxdegree', mmax_tayl);

            fprintf('Processing expft...\n');
            [expft, s_expft(k), m_expft(k), c_expft(:, k)] =...
                expm_mp(mp(A),...
                        'epsilon', eps('mp'),...
                        'algorithm', 'transfree',...
                        'abserr', false,...
                        'approx', 'taylor',...
                        'maxscaling', smax,...
                        'maxdegree', mmax_tayl);

            fprintf('Processing expp...\n');
            [expp, s_expp(k), m_expp(k), c_expp(:, k)] =...
                expm_mp(mp(A),...
                        'epsilon', eps('mp'),...
                        'algorithm', 'complexschur',...
                        'abserr', false,...
                        'approx', 'diagonalcheap',...
                        'maxscaling', smax,...
                        'maxdegree', mmax_diag);

            fprintf('Processing expfp...\n');
            [expfp, s_expfp(k), m_expfp(k), c_expfp(:, k)] =...
                expm_mp(mp(A),...
                        'epsilon', eps('mp'),...
                        'algorithm', 'transfree',...
                        'abserr', false,...
                        'approx', 'diagonalcheap',...
                        'maxscaling', smax,...
                        'maxdegree', mmax_diag);

            if (j == 1)
                old_d = mp.Digits();
                mp.Digits(20);
                fprintf('Processing conditioning...\n');
                condest1(k) = funm_condest1(double(A), @expm);
                mp.Digits(old_d);
            end
            condest1u(j, k) = condest1(k) * epsilon;

            mp.Digits(pp);
            fwderr_expmct(j, k) = norm(expmct - exact, 1) / norm_exact;
            fwderr_expotf(j, k) = norm(expotf - exact, 1) / norm_exact;
            fwderr_expt(j, k) = norm(expt - exact, 1) / norm_exact;
            fwderr_expft(j, k) = norm(expft - exact, 1) / norm_exact;
            fwderr_expp(j, k) = norm(expp - exact, 1) / norm_exact;
            fwderr_expfp(j, k) = norm(expfp - exact, 1) / norm_exact;

        end
    end

    save('test_expm_data_tmp.mat',...
         'hermitian_select', 'condest1u',...
         'fwderr_expmct', 'fwderr_expotf',...
         'fwderr_expt', 'fwderr_expft',...
         'fwderr_expp', 'fwderr_expfp');

    fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60);
else
    load test_expm_data_tmp.mat
end

%% Plot results
if plot_out
    print_legend = false;
    for j = 1:P

        Tcolors = [color_expmct; color_expotf;...
                   color_expt; color_expft; color_expp; color_expfp];
        Tstyles = {ls_expmct; ls_expotf;...
                   ls_expt; ls_expft; ls_expp; ls_expfp};
        Tmarkers = {marker_expmct; marker_expotf;...
                    marker_expt; marker_expft; marker_expp; marker_expfp};

        legend_perm = [1 2 3 4 5 6];

        indices1 = find(~hermitian_select);
        [~, perm1] = sort(condest1u(j,indices1), 'descend');
        indices2 = find(hermitian_select);
        [~, perm2] = sort(condest1u(j,indices2), 'descend');

        % Plot performance profile
        prefix_string = '';
        produce_perfprof(fwderr_expmct(j,indices1), fwderr_expotf(j,indices1),...
                         fwderr_expt(j,indices1), fwderr_expft(j,indices1),...
                         fwderr_expp(j,indices1), fwderr_expfp(j,indices1),...
                         p_vec(j), Tcolors, Tstyles, print_legend, legend_perm,...
                         png_out, tikz_out,...
                         sprintf('%snoherm_multi_%04d', prefix_string, p_vec(j)));

        produce_error_plot(condest1u(j,indices1), perm1, ...
                           fwderr_expmct(j,indices1), fwderr_expotf(j,indices1),...
                           fwderr_expt(j,indices1), fwderr_expft(j,indices1),...
                           fwderr_expp(j,indices1), fwderr_expfp(j,indices1),...
                           p_vec(j), ymins(j), ymaxs(j),...
                           Tcolors, color_cond, Tmarkers, ls_cond, msize, lw,...
                           true, print_legend, legend_perm,...
                           'noherm', png_out, tikz_out,...
                           sprintf('%snoherm_multi_%04d', prefix_string, p_vec(j)));

        produce_error_plot(condest1u(j,indices2), perm2, ...
                           fwderr_expmct(j,indices2), fwderr_expotf(j,indices2),...
                           fwderr_expt(j,indices2), fwderr_expft(j,indices2),...
                           fwderr_expp(j,indices2), fwderr_expfp(j,indices2),...
                           p_vec(j), ymins(j), ymaxs(j),...
                           Tcolors, color_cond, Tmarkers, ls_cond, msize, lw,...
                           0, false, false,...
                           'herm', png_out, tikz_out,...
                           sprintf('%sherm_multi_%04d', prefix_string, p_vec(j)));

        if (paused)
            pause;
        end

    end
end