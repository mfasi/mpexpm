function produce_error_plot(condest1u, perm, fwderr_expmct, fwderr_expotf,...
                            fwderr_expt, fwderr_expft,...
                            fwderr_expp, fwderr_expfp,...
                            prec, ymin, ymax, ...
                            Tcolors, color_cond, Tmarkers,...
                            ls_cond, msize, lw,...
                            print_yticks, print_legend, legend_perm,...
                            xlabel_string, png_out, tikz_out, filename_suffix)

    fwderr_expmct(fwderr_expmct > ymax) = ymax;
    fwderr_expotf(fwderr_expotf > ymax) = ymax;
    fwderr_expt(fwderr_expt > ymax) = ymax;
    fwderr_expft(fwderr_expft > ymax) = ymax;
    fwderr_expp(fwderr_expp > ymax) = ymax;
    fwderr_expfp(fwderr_expfp > ymax) = ymax;

    fwderr_expmct(fwderr_expmct < ymin) = ymin;
    fwderr_expmct(fwderr_expotf < ymin) = ymin;
    fwderr_expt(fwderr_expt < ymin) = ymin;
    fwderr_expft(fwderr_expft < ymin) = ymin;
    fwderr_expp(fwderr_expp < ymin) = ymin;
    fwderr_expfp(fwderr_expfp < ymin) = ymin;

    % forward error
    condest1u(isnan(condest1u)) = Inf;
    figure
    clf;
    mats = 1:1:length(condest1u);
    hold on;
    mp_semilogy(mats, condest1u(perm),...
                [], ls_cond, 'Color', color_cond, 'MarkerSize', msize,...
                'Linewidth', lw);
    mp_semilogy(mats, fwderr_expmct(perm),...
                [], 0, Tmarkers{1}, 'Color', Tcolors(1,:), 'MarkerSize', msize);
    mp_semilogy(mats, fwderr_expotf(perm),...
                [], 0, Tmarkers{2}, 'Color', Tcolors(2,:), 'MarkerSize', msize);
    mp_semilogy(mats, fwderr_expt(perm),...
                [], 0, Tmarkers{3}, 'Color', Tcolors(3,:), 'MarkerSize', msize);
    mp_semilogy(mats, fwderr_expft(perm),...
                [], 0, Tmarkers{4}, 'Color', Tcolors(4,:), 'MarkerSize', msize);
    mp_semilogy(mats, fwderr_expp(perm),...
                [], 0, Tmarkers{5}, 'Color', Tcolors(5,:), 'MarkerSize', msize);
    mp_semilogy(mats, fwderr_expfp(perm),...
                [0, length(condest1u) + 1, ymin, ymax], 4, Tmarkers{6},...
                'Color', Tcolors(6,:), 'MarkerSize', msize);
    set(gca,'XTick',0)
    if ~print_yticks
        set(gca,'YTickLabels',{})
    end
    mp.Digits(prec);

    if(print_legend)
        local_perm = [1, legend_perm+1];
        ghH = get(gca, 'Children');
        ghH = flipud(ghH);
        set(gca, 'Children', flipud(ghH(local_perm)));
        methods_names = {'kexpAu', 'expmct', 'expotf',...
                         'expt', 'expft', 'expp', 'expfp'};
        legend(methods_names{local_perm},...
               'Location', 'NE');
    end
    xlabel(xlabel_string);
    box on

    % save .png
    if png_out
        filename = sprintf('pngfigs/expm_mp_accuracy_%s', filename_suffix);
        saveas(gcf, filename, 'png');
    end

    % save .tikz
    if tikz_out

        filename = sprintf('figs/expm_mp_accuracy_%s.tikz', filename_suffix);
        matlab2tikz(filename, 'showInfo', false,...
                    'extraTikzpictureOptions',...
                    {'trim axis left', 'trim axis right'});
    end

end