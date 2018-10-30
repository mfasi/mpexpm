function [varargout] = expm_mp(A, varargin)
%EXPM_MP  Multiprecision algorithm for the matrix exponential.
%
%   This code requires the Advanpix Multiprecision Computing Toolbox (see
%   www.advanpix.com).
%
%   [X, S, M] = expm_mp(A) computes the matrix exponential X of the square
%   matrix A.
%
%   The output parameters S and M are the number of matrix products computed
%   during the squaring phase and the number of matrix multiplications (and
%   matrix inversions) required to evaluate the Pade approximant.
%
%   [...] = expm_mp(...,'precision',DIGITS) specifies the number of digits to be
%   used in the computation. Default is mp.Digits() if A is of class mp. The
%   computation is performed in single or double arithmetic if A is of class
%   single or double, respectively, and the precision is not specified.
%
%   [...] = expm_mp(...,'epsilon',EPSILON) specifies the tolerance to be used to
%   evaluate the Pade approximants. Default is machine epsilon of the precision
%   of A if A is of class 'single' or 'double', mp.eps() if A is of class 'mp'.
%
%   [...] = expm_mp(...,'maxscaling',MAXSQRTM) specifies the maximum number of
%   matrix products allowed during the squaring phase. Default is 100.
%
%   [...] = expm_mp(...,'maxdegree',MAXDEGREE) specifies the maximum degree of
%   the Pade approximants. Default is 100.
%
%   [...] = expm_mp(...,'algorithm', KIND), specifies whether the algorithm
%   should be performed on the full matrix (KIND='transfree'), on the upper
%   triangular factor of the complex Schur decomposition of the matrix (KIND =
%   'complexschur'), or on the upper quasi-triangular factor of the real Schur
%   decomposition of the matrix. If the input matrix has complex entries,
%   'complexschur' is used instead of 'realschur'. Default is 'transfree'.
%
%   [...] = expm_mp(...,'approx', KIND), where KIND='diagonal' or KIND='taylor',
%   specifies the kind of Pade approximant to use. Default is 'diagonal'. The
%   efficient evaluation of truncated Taylor series requires storing several
%   additional matrices having the same size as A.
%
%   Reference:
%   Algorithm 4.1 in M. Fasi and N. J. Higham, An Arbitrary Precision Scaling
%   and Squaring Algorithm for the Matrix Exponentia. Technical Report 201X.XX,
%   Manchester Institute for Mathematical Sciences, The University of
%   Manchester, UK, Feb XXXX.

% Copyright (c) 2017-2018, Massimiliano Fasi
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%   * Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer.
%   * Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
% EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% Check whether the Advanpix MCT toolbox is available.
% if ~check_tb_available('Advanpix Multiprecision Computing Toolbox')
%   err_msg = ['The function EXPM_MP requires the Advanpix '...
%     'Multiprecision Computing Toolbox,\nwhich does not appear '...
%     'to be installed on this system.\nSee www.advanpix.com.'];
%   error('logm_mp:missing_toolbox', err_msg);
% end

%% Parse and validate input.
    p = inputParser;
    addRequired(p, 'A', @(x)(ismatrix(x) && isnumeric(x)));
    addParameter(p, 'precision', [], @(x)(x == round(x) && x > 0))
    addParameter(p, 'epsilon', [], @(x)(x > 0 && x < 1));
    addParameter(p, 'maxscaling', 100, @(x)(x == round(x) && x > 0));
    addParameter(p, 'maxdegree', 500, @(x)(x == round(x) && x > 0));
    addParameter(p, 'algorithm', 'transfree',...
                 @(x)(ischar(x) && strcmp(x, 'transfree')...
                      || strcmp(x, 'realschur') || strcmp(x, 'complexschur')));
    addParameter(p, 'approx', 'taylor',...
                 @(x)(ischar(x) && strcmp(x, 'taylor') ||...
                      strcmp(x, 'diagonal') ||...
                      strcmp(x, 'diagonalcheap')));
    addParameter(p, 'abserr', false, @(x)(x == true || x == false));
    addParameter(p, 'timing', false, @(x)(x == true || x == false));

    parse(p,A,varargin{:});

    A = p.Results.A;
    digits = p.Results.precision;
    epsilon = p.Results.epsilon;
    maxscaling = p.Results.maxscaling;
    maxdegree = p.Results.maxdegree;
    abserr = p.Results.abserr;
    timing = p.Results.timing;

    if strcmp(p.Results.approx, 'taylor')
        usetaylor = true;
        opt_deg_pade = @opt_degrees_taylor;
        eval_pade = @expm_taylor;
        eval_pade_error = @scalar_error_taylor;
        alpha_pade = @(s, m)(alpha(s, m, 0));
    else
        usetaylor = false;
        if strcmp(p.Results.approx, 'diagonalcheap')
            opt_deg_pade = @opt_degrees_diagonalcheap;
        else
            opt_deg_pade = @opt_degrees_diagonal;
        end
        eval_pade = @expm_diagonal;
        eval_pade_error = @scalar_error_diagonal;
        alpha_pade = @(s, m)(alpha(s, m, m));
    end

    alphas = zeros(1,maxdegree);
    alphas(1) = norm(A,1);

    [n2, n] = size(A);
    if n ~= n2
        error('The matrix ``A'' must be square.');
    end

    %% Parse and validate output.
    if nargout > 4
        error('This function returns at most four values.');
    end

    %% Determine whether the computation will be single, double or mp.
    current_digits = mp.Digits();
    if ~isempty(digits) % working precision specified
        mp.Digits(digits)
        A = mp(A);
    else
        if isa(A, 'double')
            digits = 16;
        elseif isa(A, 'single')
            digits = 8;
        else
            digits = mp.Digits();
        end
        mp.Digits(digits);
    end
    if isempty(epsilon)
        epsilon = myeps(class(A));
    end
    %     epsilon = epsilon;
    curr_epsilon = epsilon;

    %% Form complex Schur form (if required).
    compute_schur = false;
    schur_string = 'complex';
    recompute_diag_blocks = true;
    switch p.Results.algorithm
      case 'transfree'
        X = A;
        if ishermitian(A)
            [V, D] = eig(A);
            if nargout > 0
                [varargout{1}] = V * diag(exp(diag(D))) / V;
            end
            if nargout > 1
                [varargout{2}] = 0;
            end
            if nargout > 2
                [varargout{3}] = 0;
            end
            if nargout > 3
                [varargout{4}] = 0;
            end
            return
        end
        if matlab.internal.math.isschur(A)
            recompute_diag_blocks = true;
        else
            recompute_diag_blocks = false;
        end
      case 'realschur'
        if matlab.internal.math.isschur(A)
            X = A;
        else
            compute_schur = true;
            schur_string = 'real';
        end
      case 'complexschur'
        if istriu(A)
            X = A;
        else
            compute_schur = true;
            schur_string = 'complex';
        end
    end

    % time(1): Computing the Schur decomposition.
    % time(2): Evaluating the bound.
    % time(3): Evalutating the Pade approximant.
    % time(4): Squaring phase.
    time = zeros(1,4);

    if compute_schur
        schur_time = tic;
        [Q, X] = schur(A, schur_string);
        time(1) = toc(schur_time);
    end
    %   if ishermitian(A)
    %     [Q, X] = eig(A);
    %   end

    % Shift the input matrix -- if needed.
    % useshift = false;
    mu = trace(X) / n;
    if abs(mu) > 10
        useshift = true;
        X(1:n+1:n^2) = X(1:n+1:n^2) - mu;
        if real(mu) >= 0
            positive_shift = true;
        else
            positive_shift = false;
        end
    else
        useshift = false;
        mu = 0;
    end

    counter = zeros(2,1);
    if isdiag(X)       % Check if T is diagonal.
        Y = diag(exp(diag(X)));
        s = 0;
        m = 0;
        if usetaylor
            currcost = 2;
        else
            currcost = 1;
        end
        degrees = [0];
    else

        X0 = X;

        % Create the vector of optimal degrees for current approximant.
        % Update maxdegree appropriately.
        degrees = opt_deg_pade(maxdegree);
        maxdegree = degrees(end);        % max degree
        currcost = 3;                    % current cost
        m = degrees(currcost);           % current degree

        %% Try all degrees lower than maxdegree.
        s = 0;

        eval_time = tic;
        if usetaylor
            factorials = factorial(mp(0:maxdegree));
            factorials_double = double(factorials);
            Xsq_powers = zeros(n,n,ceil(sqrt(maxdegree))+1,class(X0));
            Xsq_powers_double = zeros(n,n,floor(sqrt(maxdegree)),'double');
            Xsq_powers(1:n+1:n^2) = 1;
            Xsq_powers(:,:,2) = X0;
            Xsq_powers_double(:,:,1:2) = double(Xsq_powers(:,:,1:2));
        else
            Xsq_powers = zeros(n,n,ceil(sqrt(2*maxdegree))+1,class(X0));
            Xsq_powers(1:n+1:n^2) = 1;
            Xsq_powers(:,:,2) = X0 * X0;
        end
        Xsq_powers_ll = 2;
        time(3) = toc(eval_time);

        found_degree = false;
        p_factor = 1.1;
        bound_time = tic;
        tempnormexpm1 = Inf;
        compute_normexpm = true;
        cost_step = 1;

        if norm(A,1) > 1e7
            ext_prec = false;
            if usetaylor
                ext_prec = false;
                [a, Xsq_powers, tempcondq, tempnormexpm1] =...
                    eval_pade_error(Xsq_powers, alpha_pade(s, maxdegree), s, maxdegree, ext_prec);
                while (a > epsilon * tempnormexpm1 || ~isfinite(tempnormexpm1)) && s <= maxscaling
                    s = s + 1;
                    X = X / 2;
                    compute_normexpm = true;
                    [a, Xsq_powers, tempcondq, tempnormexpm1] = eval_pade_error(Xsq_powers, alpha_pade(s, maxdegree), s, maxdegree, ext_prec);
                end
            end
        end

        ext_prec = true;
        [a, Xsq_powers, tempcondq, tempnormexpm1] =...
            eval_pade_error(Xsq_powers, alpha_pade(s, m), s, m, ext_prec);
        time(2) = toc(bound_time);
        if ~abserr
            curr_epsilon = epsilon * tempnormexpm1;
        end
        old_tempcondq = tempcondq;
        bound_time = tic;
        while(~isfinite(tempnormexpm1))
            disp(tempnormexpm1)
            currcost = currcost + cost_step;
            s = s + 1;
            compute_normexpm = true;
            m = degrees(currcost);
            [a, Xsq_powers, tempcondq, tempnormexpm1] = eval_pade_error(Xsq_powers, alpha_pade(s, m), s, m, ext_prec);
        end
        time(2) = time(2) + toc(bound_time);
        a_old = Inf;
        tempnormexpm1_old = 1;
        while ~found_degree && m < maxdegree && s < maxscaling
            normqinv_bound = epsilon^(-1/8);
            if a <= curr_epsilon && tempcondq < normqinv_bound
                found_degree = true;
            elseif tempcondq >= normqinv_bound ||...
                    (a_old > 1 && abs(a_old) < abs(a)^2)
                assert(~(tempcondq >= normqinv_bound && usetaylor))
                X = X / 2;
                s = s + 1;
                compute_normexpm = true;
            else
                currcost = currcost + cost_step;
                m = degrees(currcost);
            end
            a_old = a;
            tempnormexpm1_old = tempnormexpm1;
            bound_time = tic;
            [a, Xsq_powers, tempcondq, tempnormexpm1] = eval_pade_error(Xsq_powers, alpha_pade(s, m), s, m, ext_prec);
            if ~abserr
                curr_epsilon = epsilon * tempnormexpm1;
            end
            time(2) = time(2) + toc(bound_time);
        end

        % Degree and mxm_ps are now fixed and we keep scaling until the bound
        % on the forward error is smaller than the unit roundoff.
        bound_time = tic;
        if ~found_degree
            [a, Xsq_powers, tempcondq, tempnormexpm1] = eval_pade_error(Xsq_powers, alpha_pade(s, maxdegree), s, maxdegree, ext_prec);
            if ~abserr
                curr_epsilon = epsilon * tempnormexpm1;
            end
            %         while (~isfinite(a) || a >= curr_epsilon || tempcondq >= normqinv_bound) && s < maxscaling
            while (~isfinite(a) || a >= curr_epsilon) && s < maxscaling
                X = X / 2;
                s = s + 1;
                compute_normexpm = true;
                [a, Xsq_powers, tempcondq, tempnormexpm1] = eval_pade_error(Xsq_powers, alpha_pade(s, maxdegree), s, m, ext_prec);
                if ~abserr
                    curr_epsilon = epsilon * tempnormexpm1;
                end
            end
        end
        time(2) = time(2) + toc(bound_time);

        %% Compute matrix exponential.
        eval_time = tic;
        X = X0/(2^s);
        Y = eval_pade(X, Xsq_powers, s, m);
        if recompute_diag_blocks
            Y = recompute_diagonals(X, Y);
        end
        if useshift && ~positive_shift
            Y = exp(2^-s*mu) * Y;
            X(1:n+1:n^2) = X(1:n+1:n^2) + 2^-s*mu;
        end
        time(3) = time(3) + toc(eval_time);

        %% Squaring.
        % Squaring phase.
        squaring_time = tic;
        for t = 1:s
            Y = Y * Y;
            if recompute_diag_blocks
                X = 2 * X;
                Y = recompute_diagonals(X, Y);
            end
        end
        time(4) = toc(squaring_time);

    end

    if compute_schur
        Y = Q * Y * Q';
    end
    if isreal(A)
        Y = real(Y);
    end

    if useshift && positive_shift
        Y = exp(mu) * Y;
    end

    %% Prepare output.
    if nargout > 0
        [varargout{1}] = Y;
    end
    if nargout > 1
        [varargout{2}] = s;
    end
    if nargout > 2
        [varargout{3}] = m;
    end
    if nargout > 3
        %         [varargout{4}] = counter;
        [varargout{4}] = time;
    end

    %% Restore mp working precision.
    mp.Digits(current_digits);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                             SUBFUNCTIONS                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Scalar bound evaluation and Pade approximation.

    function [e, Xsq_powers, tempcondq, tempnormexpm] = scalar_error_taylor(Xsq_powers, x, s, m, ext_prec)
    %SCALAR_ERROR_TAYLOR   Error in Taylor approx. to the exponential.
    %   SCALAR_ERROR_TAYLOR(A,X,M) computes an upper bound on the
    %   truncation error of the Taylor approximant of degree M to EXP(A).

        time_eval = tic;
        ss = ceil(sqrt(m));
        for i = Xsq_powers_ll+1:ss+1
            %             assert(all(all(Xsq_powers(:,:,i) == 0)))
            Xsq_powers(:,:,i) = Xsq_powers(:,:,i-1) * Xsq_powers(:,:,2);
        end
        Xsq_powers_double(:,:,Xsq_powers_ll+1:ss+1) = double(Xsq_powers(:,:,Xsq_powers_ll+1:ss+1));
        Xsq_powers_ll = ss+1;
        time_tmp = toc(time_eval);
        time(2) = time(2) - time_tmp;
        time(3) = time(3) + time_tmp;

        digits_old = mp.Digits();
        if ext_prec
            mp.Digits(p_factor*digits_old);
        end
        x = mp(x);
        y = sum(x.^(0:m)./factorials(1:m+1));
        e = abs(y - exp(x));
        mp.Digits(digits_old);
        tempcondq = 1;
        if compute_normexpm
            reshaped_vector = reshape(2.^(-s*(0:Xsq_powers_ll-1))./factorials_double(1:Xsq_powers_ll),[1,1,Xsq_powers_ll]);
            approx = sum(bsxfun(@times,...
                                Xsq_powers_double(:,:,1:Xsq_powers_ll),... % reshaped_matrix,...
                                reshaped_vector),3);
            tempnormexpm = norm(approx, 1);
            if abs(tempnormexpm1 - tempnormexpm) / abs(tempnormexpm) < sqrt(eps)
                compute_normexpm = false;
                %             disp('.');
            end
        else
            tempnormexpm = tempnormexpm1;
        end
    end

    function [e, Xsq_powers, tempcondq, tempnormexpm] = scalar_error_diagonal(Xsq_powers, x, s, m, ext_prec)
    %SCALAR_ERROR_DIAGONAL   Error in Pade approx. to the exponential.
    %   SCALAR_ERROR_DIAGONAL(X,M) computes an upper bound on the
    %   truncation error of the [M/M] Pade approximant to EXP(A).
        digits_old = mp.Digits();
        if ext_prec
            mp.Digits(p_factor*digits_old);
        end
        m = mp(m);
        x = mp(x);
        l_num = mp(0:1:m);
        c_num = (factorial(m)/factorial(2 * m)) ./...
                (factorial(m-l_num).*factorial(l_num)) .* factorial(2*m-l_num);
        c_den = (-1).^l_num .* c_num;
        xx = cumprod([1, repmat(x, 1, m)]);

        [Ue, Xsq_powers] = polyvalm_ps(Xsq_powers, s, double(c_num(1:2:end)), 'double');
        if m >= 1
            [Uo, Xsq_powers] =...
                polyvalm_ps(Xsq_powers, s, double(c_num(2:2:end)), 'double');
            Uo = double(X0) / 2^s * Uo;
        else
            Uo = zeros(size(X0), 'like', Ue);
        end

        tempq = Ue - Uo;

        tempqinv = inv(double(tempq));
        tempnormqinv = norm(tempqinv, 1);
        tempcondq = norm(tempq, 1) * tempnormqinv;
        tempnormexpm = norm(2 * (tempq \ Uo) + eye(size(X0), 'like', tempq));

        e = tempnormqinv * abs(exp(x) * sum(c_den .* xx) - sum(c_num .* xx));
        mp.Digits(digits_old);

    end

    function [c,mv] = mynormAm(X,p)
    %NORMAM   Estimate of 1-norm of power of matrix.
    %   NORMAM(A,m) estimates norm(A^m,1).
    %   If A has nonnegative elements the estimate is exact.
    %   [C,MV] = NORMAM(A,m) returns the estimate C and the number MV of
    %   matrix-vector products computed involving A or A^*.

    %   Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
    %   Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
    %   970-989, 2009.

        tt = 1;
        n1 = length(X);
        if isequal(X,abs(X))
            e = ones(n1,1);
            for j=1:p         % for positive matrices only
                e = X'*e;
            end
            c = norm(e,inf);
            mv = p;
        else
            % Use powers of X in Xsq_powers in order to efficienlty compute
            % A^m * y.
            if usetaylor
                mult = zeros(1, Xsq_powers_ll);
                p_dec = p;
                ll_dec = min(Xsq_powers_ll, p+1);
                while p_dec > 0
                    mult(ll_dec) = floor(p_dec / (ll_dec-1));
                    p_dec = mod(p_dec,(ll_dec-1));
                    ll_dec = min(ll_dec - 1, p_dec + 1);
                end
            else
                ll = size(Xsq_powers, 3);
                mult = zeros(1, ll);
                p_dec = p;
                ll_dec = ll;
                while p_dec > 1 && ll_dec > 1
                    mult(ll_dec) = floor(p_dec / (2*(ll_dec-1)));
                    %     p_dec = p_dec - mult(ll_dec) * 2 * (ll_dec-1)
                    p_dec = mod(p_dec,2*(ll_dec-1));
                    ll_dec = min(ll_dec - 1, floor(p_dec/2) + 1);
                end
                mult(1) = p_dec;
            end

            [c,~,~,it] = normest1(@afun_power,tt);
            mv = it(2)*tt*p;
        end

        function Z = afun_power(flag,Y)
        %AFUN_POWER  Function to evaluate matrix products needed by NORMEST1.

            if isequal(flag,'dim')
                Z = n1;
            elseif isequal(flag,'real')
                Z = isreal(A);
            else
                [n2,~] = size(Y);
                if n2 ~= n1
                    error('Dimension mismatch')
                end
                if isequal(flag,'notransp')
                    for i = find(mult(2:end))
                        for ii = 1:mult(i+1)
                            Y = double(Xsq_powers(:,:,i+1)) * Y;
                        end
                    end
                    if mult(1) ~= 0
                        for ii = 1:mult(1)
                            Y = X * Y;
                        end
                    end
                elseif isequal(flag,'transp')
                    for i = find(mult(2:end))
                        for ii = 1:mult(i+1)
                            Y = double(Xsq_powers(:,:,i+1)') * Y;
                        end
                    end
                    if mult(1) ~= 0
                        for ii = 1:mult(1)
                            Y = X' * Y;
                        end
                    end
                end

                Z = Y;

            end

        end
    end

    function y = alpha(s, k, m)
    %ALPHA    Estimate max(||X^p||^(1/p), ||X^(p+1)||^(1/(p+1))).
    %   ALPHA(A,K,M) estimates max(||X^p||^(1/p), ||X^(p+1)||^(1/(p+1)))
    %   where p is the largest integer such that p(p-1) < K+M+1.
        j = floor((1 + sqrt(4 * (m + k) + 5)) / 2);

        % Experimental strategy
        if alphas(j+1) == 0 % The second term is not known, we try the bound.
            if alphas(j) == 0
                alphas(j) = mynormAm(double(X0), j)^(1/j);
            end
            known = find(alphas ~= 0);
            low = min(known);
            high = max(known);
            bin_counter = 0;
            found_upper_bound = false;
            while low < high
                if (low + high == j+1)
                    if (alphas(j) > alphas(low) * alphas(high))
                        found_upper_bound = true;
                        break;
                    end
                end
                if bin_counter
                    low = low + 1;
                else
                    high = high - 1;
                end
                bin_counter = mod(bin_counter + 1, 2);
            end
            if found_upper_bound
                %             fprintf('.\n');
                y = alphas(j) / 2^s;
                %             i = j;
                return
            else
                assert(alphas(j+1) == 0)
                alphas(j+1) = mynormAm(double(X0), j+1)^(1/(j+1));
            end
        end
        [y,i] = max(alphas(j:j+1));
        y = y / 2^s;
    end

    function Y = polyvalm_tay_exp(A, s, m)
    %POLYVALM_PS    Paterson-Stockmeyer method for matrix polynomials.
    %   POLYVALM_PS(A,C) evaluates the matrix polynomial
    %       Y = C(1)*EYE(N) + C(2)*A + ...  + C(M)*A^(M-1) + C(M+1)*A^M
    %   where A is an N-by-N matrix and C is a vector of length M+1
    %   using Paterson-Stockmeyer algorithm.

        n = size(A,1);
        I = eye(size(A,1), class(A));

        ss = ceil(sqrt(m));
        if ss == 0
            Y = c(1) * I;
            return
        end
        rr = floor(m/ss);

        scaling = 2^s;

        Y = Xsq_powers(:,:,ss+1) / (scaling*m);
        for j = ss - 1:-1:1
            Y = (Y + Xsq_powers(:,:,j+1)) / (scaling*(m-ss+j));
        end

        Y = Xsq_powers(:,:,ss+1) * Y;
        for k = rr - 2:-1:1
            for j = ss:-1:1
                Y = (Y + Xsq_powers(:,:,j+1)) / (scaling * (k*ss+j));
            end
            Y = Xsq_powers(:,:,ss+1) * Y;
        end

        for j = ss:-1:1
            Y = (Y + Xsq_powers(:,:,j+1)) / (scaling * j);
        end
        Y(1:n + 1:n^2) = Y(1:n + 1:n^2) + 1;
    end

    function [Y, Xsq_powers] = polyvalm_ps(Xsq_powers, s, c, outputclass)
    %POLYVALM_PS    Paterson-Stockmeyer method for matrix polynomials.
    %   POLYVALM_PS(A,C) evaluates the matrix polynomial
    %       Y = C(1)*EYE(N) + C(2)*A + ...  + C(M)*A^(M-1) + C(M+1)*A^M
    %   where A is an N-by-N matrix and C is a vector of length M+1
    %   using Paterson-Stockmeyer algorithm.

        if nargin == 3
            outputclass = class(X0);
        end
        mm = length(c) - 1;
        I = Xsq_powers(:,:,1);

        ss = ceil(sqrt(mm));
        if ss == 0
            Y = c(1) * Xsq_powers(:,:,1);
            return
        end
        rr = floor(mm/ss);

        % Compute first ss+1 powers.
        time_eval = tic;
        for i=Xsq_powers_ll+1:ss+1
            if mod(i,2) == 1 % even power
                i_half = (i - 1) / 2 + 1;
                Xsq_powers(:,:,i) = Xsq_powers(:,:,i_half)^2;
            else
                assert(all(all(Xsq_powers(:,:,i) == zeros(n,n))));
                Xsq_powers(:,:,i) = Xsq_powers(:,:,2) * Xsq_powers(:,:,i-1);
            end
        end
        Xsq_powers_ll = max(Xsq_powers_ll,ss+1);
        time_tmp = toc(time_eval);
        time(2) = time(2) - time_tmp;
        time(3) = time(3) + time_tmp;

        mpowers = Xsq_powers;
        for i = 1:ss+1
            mpowers(:,:,i) = mpowers(:,:,i)/(2^(2*s*(i-1)));
        end
        mpowers = cast(mpowers, outputclass);

        % Evaluate last polynomial of degree m mod ss.
        digits_old = mp.Digits();
        mp.Digits(1.2 * digits_old);
        B = mp(c(mm+1) * mpowers(:,:,mm-ss*rr+1));
        for j=mm-1:-1:ss*rr
            if j == ss*rr
                B = B + c(ss*rr+1)*I;
            else
                B = B + c(j+1) * mpowers(:,:,mm-ss*rr-(mm-j)+1);
            end
        end
        mp.Digits(digits_old);

        % Evaluate polynomials of degree ss-1 and evaluate main polynomial using
        % Horner's method.
        Y = cast(B, outputclass);
        for kk=rr-1:-1:0
            mp.Digits(1.2 * digits_old);
            B = zeros(size(X0,1), outputclass);
            B = B + c(ss*kk+1) * I;
            for j=1:ss-1
                B = B + c(ss*kk+j+1) * mpowers(:,:,j+1);
            end
            mp.Digits(digits_old);
            Y = Y * mpowers(:,:,ss+1) + cast(B, outputclass);
        end

    end

    function S = expm_taylor(X, Xsq_powers, s, m)
    %LOGM_TAYLOR   Taylor approximation to matrix exponential.
    %   LOGM_TAYLOR(A,M) computes the Taylor approximant of degree M to
    %   LOG(EYE(SIZE(A))+A) using the Paterson-Stockmeyer algorithm.
    %     digits_old = mp.Digits();
    %     mp.Digits(2*digits_old);
    %     c = 1./factorials(1:m+1);
    %     mp.Digits(digits_old);
    %
        S = polyvalm_tay_exp(X, s, m);
    end

    function S = expm_diagonal(A, Xsq_powers, s, m)
    %LOGM_TAYLOR   Taylor approximation to matrix exponential.
    %   LOGM_TAYLOR(A,M) computes the Taylor approximant of degree M to
    %   LOG(EYE(SIZE(A))+A) using Paterson-Stockmeyer algorithm.
        digits_old = mp.Digits();
        mp.Digits(2*digits_old);
        m = mp(m);
        l_num = mp(0:1:m);
        c_num = (factorial(m)/factorial(2 * m)) ./...
                (factorial(m-l_num).*factorial(l_num)) .*...
                factorial(2*m - l_num);
        c_den = (-1).^l_num .* c_num;

        mp.Digits(digits_old);

        [Ue, ~] = polyvalm_ps(Xsq_powers, s, c_num(1:2:end));
        if m >= 1
            [Uo, ~] = polyvalm_ps(Xsq_powers, s, c_num(2:2:end));
            Uo = (X0 / 2^s) * Uo;
        else
            Uo = zeros(size(A), class(A));
        end

        if (strcmp(p.Results.approx, 'diagonalcheap'))
            Qm = (Ue - Uo);
            S = Qm \ (2*Uo) + eye(size(A), class(A));
        else
            Pm = polyvalm_ps(Xsq_powers, s, c_num);
            Qm = polyvalm_ps(Xsq_powers, s, c_den);
            S = Qm \ Pm;
        end

    end

    function degs = opt_degrees_taylor(nmax)
        degs = [1,    2,    4,    6,    9,   12,   16,   20,   25,...
                30,   36,   42,   49,   56,   64,   72,   81,   90,  100,...
                110,  121,  132,  144,  156,  169,  182,  196,  210,  225,...
                240,  256,  272,  289,  306,  324,  342,  361,  380,  400,...
                420,  441,  462,  484,  506,  529,  552,  576,  600,  625,...
                650,  676,  702,  729,  756,  784,  812,  841,  870,  900,...
                930,  961,  992, 1024, 1056, 1089, 1122, 1156, 1190, 1225,...
                1260, 1296, 1332, 1369, 1406, 1444, 1482, 1521, 1560, 1600,...
                1640, 1681, 1722, 1764, 1806, 1849, 1892, 1936, 1980, 2025,...
                2070, 2116, 2162, 2209, 2256, 2304, 2352, 2401, 2450, 2500];
        degs = degs(1:find(degs <= nmax, 1, 'last'));
    end

    function degs = opt_degrees_diagonal(nmax)
        degs = [1,    2,    3,    4,    6,    8,   10,   12,   15,...
                18,   21,   24,   28,   32,   36,   40,   45,   50,   55,...
                60,   66,   72,   78,   84,   91,   98,  105,  112,  120,...
                128,  136,  144,  153,  162,  171,  180,  190,  200,  210,...
                220,  231,  242,  253,  264,  276,  288,  300,  312,  325,...
                338,  351,  364,  378,  392,  406,  420,  435,  450,  465,...
                480,  496,  512,  528,  544,  561,  578,  595,  612,  630,...
                648,  666,  684,  703,  722,  741,  760,  780,  800,  820,...
                840,  861,  882,  903,  924,  946,  968,  990, 1012, 1035,...
                1058, 1081, 1104, 1128, 1152, 1176, 1200, 1225, 1250, 1275,...
                1300, 1326, 1352, 1378, 1404, 1431, 1458, 1485, 1512, 1540,...
                1568, 1596, 1624, 1653, 1682, 1711, 1740, 1770, 1800, 1830,...
                1860, 1891, 1922, 1953, 1984, 2016, 2048, 2080, 2112, 2145,...
                2178, 2211, 2244, 2278, 2312, 2346, 2380, 2415, 2450, 2485,...
                2520];
        degs = degs(1:find(degs >= nmax, 1, 'first'));
    end

    function degs = opt_degrees_diagonalcheap(nmax)
        degs = [1,    2,    3,    5,    7,    9,   13,   17,   21,...
                25,   31,   37,   43,   49,   57,   65,   73,   81,   91,...
                101,  111,  121,  133,  145,  157,  169,  183,  197,  211,...
                225,  241,  257,  273,  289,  307,  325,  343,  361,  381,...
                401,  421,  441,  463,  485,  507,  529,  553,  577,  601,...
                625,  651,  677,  703,  729,  757,  785,  813,  841,  871,...
                901,  931,  961,  993, 1025, 1057, 1089, 1123, 1157, 1191,...
                1225, 1261, 1297, 1333, 1369, 1407, 1445, 1483, 1521, 1561,...
                1601, 1641, 1681, 1723, 1765, 1807, 1849, 1893, 1937, 1981,...
                2025, 2071, 2117, 2163, 2209, 2257, 2305, 2353, 2401, 2451,...
                2501];
        degs = degs(1:find(degs >= nmax, 1, 'first'));
    end

    %% Recomputation of diagonal blocks.
    function F = recompute_diagonals(T, F)

        i = 1;
        while i <= size(T, 1)
            if (n == i + 1) || ((i <= n - 2) && (T(i+2,i+1) == 0))
                % start of 2-by-2 block
                if T(i+1,i) == 0 % triangular block
                    F(i:i+1,i:i+1) = expm2by2_tri(T(i:i+1,i:i+1));
                else % full block
                    F(i:i+1,i:i+1) = expm2by2_full(T(i:i+1,i:i+1));
                end
                i = i + 2; % skip second element of 2-by-2 block
            else % start of 1-by-1 block, no superdiagonal
                F(i,i) = exp(T(i,i));
                i = i + 1;
            end
        end
    end

    function F = expm2by2_full(T)
    %EXPM2BY2    Exponential of 2-by-2 matrix.
    %   EXPM2BY2(T) computes the exponential of the 2-by-2 matrix T.

    % Equation (2.2) of
    % A. H. Al-Mohy and N. J. Higham. A new scaling and squaring
    % algorithm for the matrix exponential. SIAM J. Matrix Anal. & App.
    % 31.3 (2009):970-989. DOI: 10.1137/09074721X.

    % Start of a 2x2 quasi-triangular (full) block.
        t11 = T(1,1);
        t12 = T(1,2);
        t21 = T(2,1);
        t22 = T(2,2);
        delta = sqrt((t11-t22)^2 + 4*t12*t21)/2;
        expad2 = exp((t11+t22)/2);
        coshdelta = cosh(delta);
        sinchdelta = sinch(delta);
        F(1,1) = expad2 .* (coshdelta + (t11-t22)./2.*sinchdelta);
        F(2,1) = expad2 .* t21 .* sinchdelta;
        F(1,2) = expad2 .* t12 .* sinchdelta;
        F(2,2) = expad2 .* (coshdelta + (t22-t11)./2.*sinchdelta);

    end

    function F = expm2by2_tri(T)
    %EXPM2BY2    Exponential of 2-by-2 upper triangular matrix.
    %   EXPM2BY2(T) computes the exponential of the 2-by-2 upper
    %   triangular matrix T.

    % Equation (10.42) of
    % N. J. Higham. Functions of Matrices. SIAM, 2008.
        dT = diag(T);
        expT = exp(dT);
        F = diag(expT);
        exp_arg = (dT(1)+dT(2))/2;
        sinch_arg = (dT(2)-dT(1))/2;

        if max(exp_arg,abs(sinch_arg)) < log(realmax) % Guard against overflow
            F(1,2) = T(1,2) * exp(exp_arg) * sinch(sinch_arg);
        else
            % Numerical cancellation if dT(2) ~ dT(1)
            F(1,2) = T(1,2) * (exp(dT(2))-exp(dT(1))) / (dT(2)-dT(1));
        end

    end

    function y = sinch(x)
    %SINCH    Compute sinh(x)/x.
    %   SINCH(X) returns 1 for X = 0 and SINH(X)/X otherwise.
        if x == 0
            y = 1;
        else
            y = sinh(x)/x;
        end
    end

    %% Multiprecision Computing Toolbox utility functions.

    function e = myeps(curr_class)
    %COMP_EPS    Machine epsilon.
    %   COMP_EPS(CLASS) computes the machine epsilon of class CLASS.
        if(strcmp(curr_class, 'mp'))
            e = mp('eps');
        else
            e = eps(curr_class);
        end
    end

    %% Miscellanea

    function isavailable = check_tb_available(tb_name)
    %CHECK_TB_AVAILABLE    Check availability of toolbox.
    %    CHECK_TB_AVAILABLE(TOOLBOX) checks whether a toolbox whose name
    %    matches exactly TOOLBOX is available.
        v = ver;
        if any(strcmp(tb_name, {v.Name}))
            isavailable = true;
        else
            isavailable = false;
        end
    end

end