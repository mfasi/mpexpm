function [A, n_mats] = mygallery(k,n)
    [~, n1] = expm_testmats();
    [~, n2] = mymatrices();
    [~, n3] = gallery_getall_expm();
    n_mats = n1 + n2 + n3;
    if nargin < 1
        A = [];
        return;
    elseif nargin < 2
        n = 10;
    end
    if k <= n1
        A = full(expm_testmats(k, n));
    elseif k <= n1 + n2
        A = full(mymatrices(k - n1,n));
    elseif k <= n1 + n2 + n3
        A = full(gallery_getall_expm(k - n1 - n2, n));
    end
end