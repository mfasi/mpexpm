function [F, m, s, Cp, b ] = exptayotf (A, prec, norm_max, s_pow2, shift)
% function [F, m, s, Cp, b] = exptayotf (A, prec, norm_max, s_pow2, shift)
%
% Inputs
%              A:  n x n matrix to exponentiate
% Optional
%           prec:  precision up to which e^A is approximated
%                  (default 2^(-53) for double or higher precision data,
%                  2^(-24) for single precision data)
%       norm_max:  norm up to which we accept A without scaling
%                  (default 3.5 for double or higher precision data,
%                  6.3 for single precision data)
%         s_pow2:  if true, the scaling is a power of 2
%                  (default false)
%          shift:  if true, shift is performed
%                  (default true)
%
% Outputs
%              F:  approximation of the matrix exponential of A
%              m:  degree of approximation used
%              s:  log2 of the employed scaling parameter
%             Cp:  total matrix by matrix products
%              b:  vector containing rho values
%

% Version update: 15 / May / 2018.
% Copyright: Marco Caliari and Franco Zivcovich. This function is
% provided as it.

% Inputs handling
allow_scal_ref = true; % by default we allow the scaling refinement
issingle = isa (A, 'single');
if ((nargin == 1) || isempty (prec))
  prec = 2 ^ (-53) * ~issingle + 2 ^ (-24) * issingle;
end
if ((nargin <= 2) || isempty (norm_max))
  norm_max = 3.5 * ~issingle + 6.3 * issingle;
end
if ((nargin <= 3) || (isempty (s_pow2)))
  s_pow2 = false; % by default we allow s to be sum of two powers of 2
end
if ((nargin <= 4) || (isempty (shift)))
  shift = true; % by default we allow the shifting technique
end
scal_ref = allow_scal_ref;
scal4deg = allow_scal_ref;

%% Entry point
%% Type handling
n = length (A);
if ((issparse (A)) && (nnz (A) > 5 * n)) % 5*n chosen empirically
  B{1} = full (A);
else
  B{1} = A;
end

%% Preprocessing Shift
b = zeros (1, 6);
b(1) = norm (B{1}, 1);
tol = b(1) * prec; % compute ||A||_1 * prec
mu = full (trace (A)) / n * shift;
if ((abs (mu) / b(1)) > 1e-3) % if mu is nonnegligible we compute ||A-mu*I||_1
  B{1}(1:n + 1:n ^ 2) = B{1}(1:n + 1:n ^ 2) - mu;
  b(1) = norm (B{1}, 1);
else % if mu is negligible we save a norm evaluation
  mu = 0;
end

%% Sudden exit
if ( b(1)^2 <= 2 * tol ) % a fast backward error analysis may grant sudden exit
  B{1}(1:n + 1:n ^ 2) = B{1}(1:n + 1:n ^ 2) + 1; % degree 1
  F = B{1} * exp (mu); m = 1; s = 0; Cp = 0;
  return
end

%% Entry Point
z = 2;
B{2} = B{1} * B{1};
b(2) = nthroot (double (norm (B{2}, 1)), 2); % double needed for symbolic input
rho = min (b(1), b(2));

%% Creation of X, n - by - 4 matrix auxiliary to 1-norm estimate algorithm
if (isnumeric (A))
  X = [ones(n, 1, class (A)), ...
       sign(rand (n, 4, class (A)) - 0.5)] / n;
else % A is vpa, e.g.
  X = [ones(n, 1), ...
       sign(rand (n, 4) - 0.5)] / n;
end

%% Degree and scaling selection
Cp = 2;
i_o = false;
s = 1;
m = 0; % m = 4 but we need m = 0 to well update X in the first while iteration
while (~i_o)
  s0 = s;
  [s, q] = scaling (rho / norm_max, s_pow2);
  
  X = B{z} * X * ( ( ( s0 / s ) ^ m ) / s^z );
  m = (Cp - z + 2) * z; % refresh m
  threshold = exp (gammaln (m + 1) + log (min (1, tol / s)));
  if ((any (isnan (X(:)))) || (any (isinf (X(:)))) || (threshold == Inf))
    scal_ref = false;
    break
  end

  %% on the fly estimate of the backward error
  i_o = bea_check (B, z, m, threshold, s, X); % i_o = (backwerr < threshold)
  Cp = Cp + ~i_o; % if threshold is exceeded we must face a greater effort
  if (ceil (Cp / 2) + 1 > z) % if this then we must increase z
    z = z + 1;
    B{z} = B{ceil (z / 2)} * B{floor (z / 2)};
    if (scal_ref)
      b(z) = nthroot (double (norm (B{z}, 1)), z);
      %% we won't try the scaling refinement
      %% if ||B^p||^(1/p) doesn't decrease fast enough
      scal_ref = (b(z) < 0.8*b(z-1)) || ...
                 ((b(z) < 0.8*b(z-2)) && (b(z) > b(z-1)));
      %% we won't trade a scaling Cp for a degree Cp
      %% if ||B^p||^(1/p) doesn't decrease fast enough
      scal4deg = (b(z) < 0.7*b(z-1)) || ...
                 ((b(z) < 0.7*b(z-2)) && (b(z) > b(z-1)));
      rho = min (rho, b(z)); % refresh rough overestimate of spectral radius
    end
  end
end


%% Scaling refinement: try to reduce s.
s0 = s;
q0 = q;
while ((s > 1) && allow_scal_ref && scal_ref && i_o )
  %% Try and reduce s
  [s, q] = scaling (rho / norm_max, s_pow2, ceil (log2 (s0 / 2)));
  X =  X * (s0 / s) ^ m;
  threshold = exp (gammaln (m + 1) + log (min (1, tol / s)));
%   threshold = tol / s;
  if ((any (isnan (X(:)))) || (any (isinf (X(:)))) || (threshold == Inf))
    break
  end
  i_o = bea_check (B, z, m, threshold, s, X); % i_o = (backwerr < threshold)
  if (i_o) % if threshold isn't exceeded we accept the maybe cheaper scaling
    [ s0, q0 ] = deal (s, q);
    continue % start now next while iteration
  end
  if (~scal4deg)
    %% if we didn't continue few lines ago means we unsuccessfully tried
    %% to reduce the scaling parameter, if we are not willing to trade
    %% a scaling Cp for a degree Cp we must quit the scaling refinement part
    break
  end
  %% Try and reduce s by increasing m
  X = B{z} * (X / s ^ z);
  threshold = exp (gammaln (m + z + 1) + log (min (1, tol / s)));
%   threshold = tol / s;
  if ((any (isnan (X(:)))) || (any (isinf (X(:)))) || (threshold == Inf))
    break
  end
  i_o = bea_check (B, z, m + z, threshold, s, X);
  Cp = Cp + i_o;
  if (ceil (Cp / 2) + 1 > z)
    z = z + 1;
    B{z} = B{ceil (z / 2)} * B{floor (z / 2)};
    if (scal_ref)
      b(z) = nthroot (double (norm (B{z}, 1)), z);
      %% we won't try the scaling refinement
      %% if ||B^p||^(1/p) doesn't decrease fast enough
      scal_ref = (b(z) < 0.8*b(z-1)) || ...
                 ((b(z) < 0.8*b(z-2)) && (b(z) > b(z-1)));
      %% we won't trade a scaling Cp for a degree Cp
      %% if ||B^p||^(1/p) doesn't decrease fast enough
      scal4deg = (b(z) < 0.7*b(z-1)) || ...
                 ((b(z) < 0.7*b(z-2)) && (b(z) > b(z-1)));
      rho = min (rho, b(z));
    end
  end
  m = (Cp - z + 2) * z;
  if (i_o) % if threshold isn't exceeded we trade a scaling Cp for a degree Cp
    [ s0, q0 ] = deal (s, q);
  end
end
[ s, q ] = deal (s0, q0);

%% Compute exp (A / s)
F = MPS (B, m, z, s, n); % Modified Paterson Stockmeyer algorithm
%% Shift back
if (shift && (real (mu) < 0)) % scale mu only if mu<0. See AMH11,\S3.1
  F = F * exp (mu / s);
end
%% Recover exp (A)
if (isinf (q))
  for i = 1:log2 (s)
    F = F * F;
    Cp = Cp + 1;
  end
else
  for i = 1:q
    F = F * F;
    Cp = Cp + 1;
  end
  A = F;
  for i = q+1:floor (log2 (s))
    F = F * F;
    Cp = Cp + 1;
  end
  F = F * A;
  Cp = Cp + 1;
end

%% Shift back
if (shift && (real(mu) >= 0)) % See AMH11,\S3.1
  F = F * exp(mu);
end
s = log2 (s); % return s according to the literature












function [s, q] = scaling (r, s_pow2, p0)
%% Determine the scaling s and the remainder q following the new scaling
%% algorithm.
if (s_pow2) % we force s to be a power of 2
  [ p, q ] = deal ( ceil (max (log2 (r), 0)), -Inf );
else
  [ p, q ] = deal ( floor (max (log2 (r), 0)), -Inf );
  if (r > 2 ^ p)
    q = max (ceil (log2 (r - 2 ^ p)), 0);
    if (p == q)
      p = p + 1;
      q = -Inf;
    end
  end
end
if ((nargin == 3) && (2 ^ p0 < 2 ^ p + 2 ^ q))
  p = p0;
  q = -Inf;
end
s = 2 ^ p + 2 ^ q;

function i_o = bea_check (B, z, m, threshold, s, X)
%% It checks whether the Paterson--Stockmeyer triangular inequality backward
%% error is smaller than the tolerance t.
hcoef = s;
est_old = 0;
for p = 0:m / z
  P = B{1} / (hcoef * (m + p * z + 1));
  hcoef = -hcoef * (p * z + 1) * s;
  for j = 2:z
    P = P + B{j} / (hcoef * (m + p * z + j));
    if (p==0 && z>2 && j<z && ( max (sum (abs (P*X))) >= threshold))
      i_o = false;
      return
    end
    hcoef = -hcoef * (p * z + j) * s;
  end
  est = lazy_est1 (@Af, B{z}, P, m / z + p, s ^ z, threshold, X);
  threshold = threshold - est;
  if (threshold <= 0)
    i_o = false;
    return
  elseif ((est <= threshold) && (est <= est_old))
    i_o = true;
    return
  end
  est_old = est;
end
%% If you're reading this something odd may be going on: reject configuration.
i_o = false; % safety reject

function F = MPS (B, m, z, s, n)
%% Modified Horner-Paterson-Stockmeyer scheme adapted to fit
%% scaled-Taylor method with z divisor of m.
F = full (B{z}) / (s * m);
for j = z - 1:-1:1
  F = (F + B {j}) / (s * (m - z + j));
end
F = B{z} * F;
for k = m / z - 2:-1:1
  for j = z:-1:1
    F = (F + B {j}) / (s * (k * z + j));
  end
  F = B{z} * F;
end
for j = z:-1:1
  F = (F + B {j}) / (s * j);
end
if ((issparse (F)) && (nnz (F) > 5 * n ))
  F = full (F);
end
F(1:n + 1:n ^ 2) = F(1:n + 1:n ^ 2) + 1;

function est = lazy_est1 (Af, B, P, pow, s, threshold, X)
%% Returns the lower estimate of 1-norm of matrix by a greedy block
%% 1-norm power method. If at any time the estimate exceeds tol then
%% it stops and returns it.
if ((nargin < 7) || any (isinf (X(:))) || any (isnan (X(:))))
  n = length (B);
  v = Af ('notransp', ones(n, 1, class(B))/n , B, P, pow, s);
  est = norm (v, 1);
else  % laziness, we save products
  v = P * X;
  est = sum (abs (v));
  [~, j] = max (est, [], 2);
  v = v(:, j);
  est = est(j);
end
if (est > threshold)
  if (strcmp ('sym', class (est)))
    % cast est to double in case it's symbolic
    est = double (est);
  end
  return
end
xi = sign (v);
x = Af ('transp', xi, B, P, pow, s);
iter = 2;
out = 0;
while ((iter < 6) && (out == 0))
  [~, j] = max (abs (x));
  v = Af ('notransp', B(:, j) / s, B, P, pow - 1, s);
  est_old = est;
  est = norm (v, 1);
  if ((est > threshold) || ...
      ((est <= est_old) || (norm (sign (v) - xi, Inf) == 0)))
    if (strcmp ('sym', class (est)))
      %% cast est to double in case it's symbolic
      est = double (est);
    end
    return % laziness
  end
  xi = sign (v);
  x = Af ('transp', xi, B, P, pow, s);
  out = ( max (abs (x)) == x(j) );
  iter = iter + 1;
end
if (strcmp ('sym', class (est)))
  est = double (est);
end

function x = Af (flag, x, A, P, i, s)
%% Function auxiliary to the norm-1 estimate.
%% Produces:
%% (A / s) ^ i * P * x (flag = 'notransp')
%% (A' / s) ^ i * P' * x (flag = 'transp')
%% length (A) (flag = 'dim')
switch flag
  case {'notransp'}
    if (~isempty (P))
      x = P * x;
    end
    for j = 1:i
      x = A * (x / s);
    end
  case {'transp'}
    if ~isempty (P),
      x = P' * x;
    end
    for j = 1:i
      x = A' * (x / s);
    end
end