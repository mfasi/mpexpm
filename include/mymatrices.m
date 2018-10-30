function [A, n_mats] = mymatrices(k,n,class)
    n_mats = 18;
    if nargin < 1
        A = [];
        return;
    elseif nargin < 2
        n = 10;
    end
    if nargin < 3
        class = 'double';
    end
    epsilon = eps();
    switch k
      case 1
        A = full(gallery('tridiag',2,-1,-2,1));
      case 2
        A = [1 1; 1 1 + 10 * epsilon];
      case 3
        A = [10 0 0; 0 1 1; 0 1 1 + 10 * epsilon];
      case 4
        A = [];
        for i = 1:4
            A = blkdiag(A,...
                        full(gallery('tridiag',i,0,i,1)));
        end
      case 5
        n = 10;
        D = diag([zeros(1, n-1) 1]);
        Q = triu(ones(n));
        A = Q * D / Q;
      case 6
        A = diag([0, 1, 1e6]);
      case 7
        A = [1e-4 0; 0 1e4];
      case 8
        A = ones(2,2);
      case 9
        n = 10;
        A = toeplitz([16-3i, (4+3i)/8, zeros(1,n-2)],[16-3i, -5, zeros(1,n-2)]);
      case 10
        A = [1 1; 0 1e2];
      case 11
        A = [1 1e3; 1e3 1];
      case 12
        A = [1 2 3; 1 2 3; 1 2 3];
      case 13
        t = -pi/2;
        A = [cos(t) -sin(t); sin(t) cos(t)];
      case 14
        v = eye(n,1);
        A = eye(n) - v*v';
      case 15
        A = [100 2 3; 4 5 6; 7 8 100];
      case 16
        A = [1 1 1; 1 1 1 + 10*eps; 1 1 1+100*eps()];
      case 17
        A = [1 2 3; 4 5 6; 7 8 1e2];
      case 18
        A = [1 1 1 0.1; 1 1 1 10*eps; 1 1 1 100*eps; 1 1 1 1000*eps];
    end
end