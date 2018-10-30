function [A,n_mats] = expm_testmats(k, n)
%EXPM_TESTMATS  Test matrices for matrix exponential.
%   [A, NMATS] = EXPM_TESTMATS(K,N) selects the K'th test matrix.
%    NMATS is the number of test matrices available.
%    N sets the dimension of those matrices for which the dimension
%    is variable.

% Authors: Awad H. Al-Mohy, Massimiliano Fasi, and Nicholas J. Higham.

    n_mats = 38;
    if nargin < 1
        A = [];
        return;
    end
    if nargin < 2,
        n = 10;
    end

    switch k
      case 1
        % \cite[Test 1]{ward77}.
        A = [4 2 0; 1 4 1; 1 1 4];
      case 2
        % \cite[Test 2]{ward77}.
        A = [29.87942128909879     .7815750847907159 -2.289519314033932
             .7815750847907159   25.72656945571064    8.680737820540137
             -2.289519314033932   8.680737820540137  34.39400925519054];
      case 3
        % \cite[Test 3]{ward77}.
        A = [-131 19 18;
             -390 56 54;
             -387 57 52];
      case 4
        % \cite[Test 4]{ward77}.
        A = gallery('forsythe',10,1e-10,0);

      case 5
        % \cite[p. 370]{naha95}.
        T = [1 10 100; 1 9 100; 1 11 99];
        A = T*[0.001 0 0; 0 1 0; 0 0 100]/T;

      case 6
        % \cite[Ex.~2]{kela98}.
        A = [0.1 1e6; 0 0.1];

      case 7
        % \cite[p.~655]{kela98}.
        A = [0  3.8e3 0    0   0
             0 -3.8e3 1    0   0
             0 0     -1  5.5e6 0
             0 0      0 -5.5e6 2.7e7
             0 0      0   0   -2.7e7];

      case 8
        % \cite[Ex.~3.10]{dipa00}
        w = 1.3; x = 1e6; n = 8;
        A = (1/n) * [w*ones(n/2) x*ones(n/2)
                     zeros(n/2)  -w*ones(n/2)];

      case 9
        A = rosser;
        A = 2.05*A/norm(A,1);  % Bad case for expm re. cost.

      case 10
        A = [0 1e4;
             -1e4 0];  % exp = [cos(x) sin(x); - sin(x) cos(x)], x = 100;

      case 11
        A = 1e2*triu(randn(n),1);  % Nilpotent.

      case 12 % log of Cholesky factor of Pascal matrix. See \cite{edst03}.
        A = zeros(n); A(n+1:n+1:n^2) = 1:n-1;

      case 13 % \cite[p.~206]{kela89}
        A = [48 -49 50 49; 0 -2 100 0; 0 -1 -2 1; -50 50 50 -52];

      case 14 % \cite[p.~7, Ex I]{pang85}
        A = [0    30 1   1  1  1
             -100   0 1   1  1  1
             0     0 0  -6  1  1
             0     0 500 0  1  1
             0     0 0   0  0  200
             0     0 0   0 -15 0];

      case 15 % \cite[p.~9, Ex II]{pang85}
              % My interpretation of their matrix for arbitrary n.
              % N = 31 corresponds to the matrix in above ref.
        A = gallery('triw',n,1);  m = (n-1)/2;
        A = A - diag(diag(A)) + diag(-m:m)*sqrt(-1);
        for i = 1:n-1, A(i,i+1) = -2*(n-1)-2 + 4*i; end

      case 16 % \cite[p.~10, Ex III]{pang85}
        A = gallery('triw',n,1,1);
        A = A - diag(diag(A)) + diag(-(n-1)/2:(n-1)/2);

      case 17
        % \cite[Ex.~5]{kela89}.
        A = [0 1e6; 0 0];   % Same as case 6 but with ei'val 0.1 -> 0.

      case 18
        % \cite[(52)]{jemc05}.
        g = [0.6 0.6 4.0]; b = [2.0 0.75];
        A = [-g(1)       0    g(1)*b(1)
             0        -g(2) g(2)*b(2)
             -g(1)*g(3)  g(3) -g(3)*(1-g(1)*b(1))];

      case 19
        % \cite[(55)]{jemc05}.
        g = [1.5 0.5 3.0 2.0 0.4 0.03]; b = [0.6 7.0];
        A1 = [-g(5)     0      0
              0      -g(1)    0
              g(4)     g(4)   -g(3)];
        A2 = [-g(6)       0    g(6)*b(2)
              0        -g(2)  g(2)*b(1)
              0         g(4) -g(4)];
        A = [zeros(3) eye(3); A2 A1];

      case 20
        % \cite[Ex.~3]{kela98}.
        A = [-1 1e7; 0 -1e7];

      case 21
        % \cite[(21)]{mopa03}.
        Thalf = [3.8235*60*24 3.10 26.8 19.9]/60;  % Half lives in seconds/
        a = log(2)./Thalf;  % decay constant
        A = diag(-a) + diag(a(1:end-1),-1);

      case 22
        % \cite[(26)]{mopa03}.
        a1 = 0.01145;
        a2 = 0.2270;
        A = [-a1              0  0
             0.3594*a1     -a2  0
             0.6406*a1     a2  0];

      case 23
        % \cite[Table 1]{kase99}.
        a = [4.916e-18
             3.329e-7
             8.983e-14
             2.852e-13
             1.373e-11
             2.098e-6
             9.850e-10
             1.601e-6
             5.796e-8
             0.000];
        A = diag(-a) + diag(a(1:end-1),-1);

      case 24
        % Jitse Niesen sent me this example.
        lambda = 1e6 * 1i;
        mu = 1/2*(-1+sqrt(1+4*lambda));
        A = [ 0, 1; lambda, -1 ] - mu*eye(2);

      case 25 % Awad

        A = [1 1e17;0 1];

      case 26 % Awad

        b = 1e3; x = 1e10;
        A = [ 1-b/2   b/2 ; -b/2   1+b/2 ];
        A = [A          x*ones(2);
             zeros(2)       -A    ];

      case 27 % Awad
        b = 1e4;
        A = [ 1-b/2   b/2 ; -b/2   1+b/2 ];

      case 28 % Awad
        b = 1e2;
        A = [ 1-b/2   b/2 ; -b^4/2   1+b/2 ];

      case 29
        % EIGTOOL.
        A = godunov_demo/100;

      case 30
        % \cite[(14.17), p. 141]{trem05}.
        A = 10*[0 1 2; -0.01 0 3; 0 0 0];

      case 31
        A = triu(schur(gallery('invol',13),'complex'),1);
      case 32
        % \cite{kuda10}
        alpha = 1; beta = 1;  % No values are given in the paper, unfortunately.
        A = -eye(n) + alpha/2*(diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));
        A(1,2) = beta; A(n,n-1) = beta;
      case 33
        % \cite[Benchmark #1]{lara17}
        % \cite[Problem 1]{zhao17}
        A = [-3.328853448977761e-07 4.915959875924379e-18;
             0    -4.915959875924379e-18];
      case 34
        % \cite[Benchmark #2]{lara17}
        % \cite[Problem 2]{zhao17}
        A = [-2.974063693062615e-07            0      1.024464026382002e-14;
             2.974063693062615e-07 -1.379680196333551e-13                 0;
             0                0     -1.024464026382002e-14];
      case 35
        % \cite[Benchmark #3]{lara17}
        % \cite[Problem 3]{zhao17}
        A = [-2.421897905520424e-03            0      5.383443102348909e-03;
             0                -3.200125487349701e-04            0;
             0                 3.200125487349701e-04 -5.398342527725431e-03];
      case 36
        % \cite[Benchmark #4]{lara17}
        % \cite[Problem 4]{zhao17}
        A = [-1.000000000000312e-04         0                  0      0;
             1.000000000000000e-04 -1.000000000009379e-04     0      0;
             0  1.000000000000000e-04 -1.188523972153541e-06  0;
             0     0          1.188523972153541e-06 -1.024464026382002e-14];
      case 37
        % \cite[Benchmark #5]{lara17}
        % \cite[Problem 5]{zhao17}
        A = sparse(...
            [ 1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6,...
              7, 7, 8, 8, 9, 9,10,10,11,11,12,12],...
            [ 1, 7, 1, 2, 2, 3,10, 3, 4, 4, 5,12, 5, 6, 6,...
              7, 7, 8, 6, 9, 8,10,10,11,11,12],...
            [ -1.000000000000049e-04
              5.880666420493406e-14
              1.000000000000000e-04
              -4.926419193745169e-04
              4.926419193745169e-04
              -3.405151448232769e-06
              2.980258985838552e-12
              3.405151448232769e-06
              -1.000000009110124e-04
              1.000000000000000e-04
              -1.000000033477380e-04
              1.212838692746004e-09
              1.000000000000000e-04
              -1.000015370544945e-04
              1.000000000000000e-04
              -1.000000000588073e-04
              1.000000000000000e-04
              -3.885005720114481e-05
              1.537023753355886e-09
              -5.077325179294990e-11
              3.885005720114481e-05
              -1.000000029802590e-04
              1.000000000000000e-04
              -1.906345381077957e-05
              1.906345381077957e-05
              -1.212838692746004e-09]);
      case 38
        % \cite[Benchmark #6]{lara17}
        % \cite[Problem 6]{zhao17}
        A = sparse([1  1  2  2  2  3  3  3  4  5  5  6  6  7  7  8  8],...
                   [1  4  1  2  4  2  3  4  4  4  5  5  6  6  7  7  8],...
                   [-2.930607054625170e-05
                    1.292290622141271e-07
                    2.446793135977101e-05
                    -2.106574217602557e-05
                    2.051479647948103e-08
                    2.106574217602557e-05
                    -9.549786402447881e-15
                    1.855074206409039e-12
                    -1.100000000000049e-04
                    1.000000000000000e-04
                    -4.926419193745169e-04
                    4.926419193745169e-04
                    -3.405151448232769e-06
                    3.405151448232769e-06
                    -1.000000091101239e-05
                    1.000000000000000e-05
                    -3.347737955438215e-12]);
        %    Easy example.
        %    case 13 % 7-by-7.  Idea from \cite[Matrix F, p. 81 ]{kags77a}
        %    e = -1; f = -1.1;
        %    A = blkdiag(compan([1, -4*e, 6*e^2, -4*e^3 e^4]), ...
        %                compan([1, -3*f, 3*f^2, -f^3]));

        %    Removed because too ill conditioned.
        %    case 8
        %    % \cite[Ex.~3]{dahi03}
        %    A = gallery('triw',4,2^(60)) - diag([17 17 2 2]);
        %    % A = A/1e4;  % Make less ill conditioned.

        %    case 8
        %    % My concoction.
        %    A = gallery('randsvd', 8, 1e14);
        %
        %    case 9
        %    % My concoction.
        %    A = randjorth(4,4,-1e8)/1e3;

    end
end