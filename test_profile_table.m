warning off

sizes = [10, 20, 50, 100, 200, 500, 1000];
nsizes = length(sizes);

p = 34;
mp.Digits(p);

filename = sprintf('tabs/table_profile_td_%04d.tex', p);
fileid_td = fopen(filename,'w');
fprintf(fileid_td, ['\\begin{tabularx}{\\textwidth}',...
    '{@{\\extracolsep{\\fill}}rr|rrrrrr|rrrrrr}\n']);
fprintf(fileid_td, '\\toprule\n');
fprintf(fileid_td, ['\\multicolumn{2}{c|}{} & ',...
    '\\multicolumn{6}{c|}{\\expft} & ',...
    '\\multicolumn{6}{c}{\\expfp} \\\\\n']);
fprintf(fileid_td, [' & $n$ & ',...
    '$M_{sqr}$ & $M_{eval}$ & $T_{bnd}$ & $T_{eval}$ & $T_{sqr}$ & $T_{tot}$ & ',...
    '$M_{sqr}$ & $M_{eval}$ & $T_{bnd}$ & $T_{eval}$ & $T_{sqr}$ & $T_{tot}$ \\\\\n']);
fprintf(fileid_td, '\\midrule\n');

mats_id = {'\verb|A|', '\verb|B|', '\verb|C|'};

nmat = 3;
for j = 1:nmat

    for i = 1:nsizes

        n = sizes(i);

        switch j
            case 1
              A = 1e2*triu(gallery('pei',n),1);
          case 2
            A = zeros(n); A(n+1:n+1:n^2) = 1:n-1;
          case 3
            A = expm(gallery('lotkin', n));
        end

        [X, s_t, m_t, time_t] = expm_mp(mp(A),...
            'approx', 'taylor',...
            'epsilon', mp('eps'),...
            'maxscaling', 100,...
            'maxdegree', 400,...
            'timing', true);

        [X, s_d, m_d, time_d] = expm_mp(mp(A),...
            'approx', 'diagonalcheap',...
            'epsilon', mp('eps'),...
            'maxscaling', 100,...
            'maxdegree', 400,...
            'timing', true);

        [~] = exptayotf(mp(A));

        s1 = sum(time_t)/100;
        s2 = sum(time_d)/100;

        if (i == 1)
            id_field = mats_id{j};
        else
            id_field = '        ';
        end
        fprintf(fileid_td,...
            '%s & %3d & %2d & %2d & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %4.1f & %2d & %2d & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %4.1f \\\\\n',...
            id_field, sizes(i),...
            sum(s_t), m_t, sum(time_t(2))/s1, time_t(3)/s1 , time_t(4)/s1,s1*100,...
            sum(s_d), m_d, sum(time_d(2))/s2, time_d(3)/s2, time_d(4)/s2,s2*100);

        fprintf('%s & %3d & %2d & %2d & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %4.1f & %2d & %2d & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %4.1f \\\\\n',...
            id_field, sizes(i),...
            sum(s_t), m_t, sum(time_t(2))/s1, time_t(3)/s1 , time_t(4)/s1,s1*100,...
            sum(s_d), m_d, sum(time_d(2))/s2, time_d(3)/s2, time_d(4)/s2,s2*100);
    end
    if (i == nsizes && j ~= nmat)
        fprintf(fileid_td, '\\midrule\n');
    end
end

fprintf(fileid_td, '\\bottomrule\n');
fprintf(fileid_td, '\\end{tabularx}');