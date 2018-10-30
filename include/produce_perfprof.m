function produce_perfprof(fwderr_expmct, fwderr_expotf,...
                          fwderr_expt, fwderr_expft, fwderr_expp, fwderr_expfp,...
                          prec, Tcolors, Tstyles, print_legend, legend_perm,...
                          png_out, tikz_out, filename_suffix)

    mp.Digits(prec);
    perfprof_norm = @(x)(max(x, x*(1-5e-2) + 5e-2*mp('eps')/2));

    % performance profiles
    figure
    clf
    T = [perfprof_norm(fwderr_expmct);...
         perfprof_norm(fwderr_expotf);...
         perfprof_norm(fwderr_expt);...
         perfprof_norm(fwderr_expft);...
         perfprof_norm(fwderr_expp);...
         perfprof_norm(fwderr_expfp)]';
    perfprof(T, 25, Tcolors, Tstyles, 1);
    axis([1 25 0 1])
    if(print_legend)
        local_perm = legend_perm;
        ghH = get(gca, 'Children');
        ghH = flipud(ghH);
        set(gca, 'Children', flipud(ghH(local_perm)));
        methods_names = {'expmct', 'expotf',...
                         'expt', 'expft', 'expp', 'expfp'};
        legend(methods_names{local_perm},...
               'Location', 'SE');
    end
    xlabel('theta');

    % save .png
    if png_out
        filename = sprintf('pngfigs/expm_mp_perfprof_%s', filename_suffix);
        saveas(gcf, filename, 'png');
    end

    % tikz
    if tikz_out
        filename = sprintf('figs/expm_mp_perfprof_%s.tikz', filename_suffix);
        matlab2tikz(filename, 'showInfo', false);
    end
end