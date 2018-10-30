function [A, n_mats] = gallery_getall_expm(k, n)

    indices = [101:130, 201:216, 301:328, 401:402];

    n_mats = length(indices);
    if nargin < 1
        A = [];
        return;
    end
    if nargin < 2,
        n = 10;
    end
    if k < 100
        local_index = indices(k);
    else
        local_index = k;
    end

    switch local_index

        % Positive real eigenvalues
      case 101, sA='cauchy';A=gallery(sA,n); %symm % optional parameters
      case 102, sA='condex';A=gallery(sA,2,4,6);
      case 103, sA='condex';A=gallery(sA,3,2);
      case 104, sA='condex';A=gallery(sA,n,3); %lower triangular
      case 105, sA='condex';A=gallery(sA,n,4,100); %symm real
      case 106, sA='dorr';A=full(gallery(sA,n,100));
      case 107, sA='dramadah';A=gallery(sA,n,2);
      case 108, sA='frank';A=gallery(sA,n);
      case 109, sA='gcdmat';A=gallery(sA,n); %symm real
      case 110, sA='grcar';A=gallery(sA,n);

      case 111, sA='hanowa';A=gallery(sA,n);
      case 112, A=hilb(n);sA='hilb'; %symm real
      case 113, sA='invhess';A=gallery(sA,n);
      case 114, sA='jordbloc';A=gallery(sA,n,1); %symm real
      case 115, sA='kahan';A=gallery(sA,n);
      case 116, sA='lehmer';A=gallery(sA,n); %symm real
      case 117, sA='minij';A=gallery(sA,n); %symm real
      case 118, sA='moler';A=gallery(sA,n); %symm real
      case 119, sA='parter';A=gallery(sA,n);
      case 120, sA='pei';A=gallery(sA,n);

      case 121, sA='poisson';A=gallery(sA,ceil(sqrt(n))); %symm real, n^2
      case 122, sA='prolate';A=gallery(sA,n,1); %symm real Toeplitz
      case 123, sA='randcorr';A=gallery(sA,n); %symm real
      case 124, sA='sampling';A=gallery(sA,n);
      case 125, sA='toeppd';A=gallery(sA,n); %symm real
      case 126, sA='tridiag';A=full(gallery(sA,n));
      case 127, sA='wathen';A=full(gallery(sA,ceil(n^(1/4)),ceil(n^(1/4)))); %symm real
      case 128, sA='wilk';A=full(gallery(sA,3));
      case 129, sA='wilk';A=full(gallery(sA,4));
      case 130, sA='wilk';A=full(gallery(sA,5)); %symm real

        % Real eigenvalues.
      case 201, sA='binomial'; A=gallery(sA,n);
      case 202, sA='fiedler';A=gallery(sA,n);
      case 203, sA = 'house'; [v, beta] = gallery(sA, n);
        A = eye(n) - beta * (v * v');
      case 204, sA='jordbloc';A=gallery(sA,n,2); %symm real
      case 205, sA='kms';A=gallery(sA,n); %symm real
      case 206, sA='lesp';A=gallery(sA,n);
      case 207, sA='lotkin';A=gallery(sA,n);
      case 208, sA='orthog';A=gallery(sA,n,1); %symm real
      case 209, sA='orthog';A=gallery(sA,n,2); %symm real
      case 210, sA='orthog';A=gallery(sA,n,5); %symm real

      case 211, sA='orthog';A=gallery(sA,n,6); %symm real
      case 212, sA='orthog';A=gallery(sA,n,-1); %symm real
      case 213, sA='redheff';A=gallery(sA,n);
      case 214, sA='riemann';A=gallery(sA,n);
      case 215, sA='ris';A=gallery(sA,n,1e1);
      case 215, sA='wilk';A=full(gallery(sA,21)); %symm real
      case 216, sA='clement';A=gallery(sA,n,1); %symm


        % Complex eigenvalues
      case 301, sA='chebspec';A=gallery(sA,n);
      case 302, sA='chebvand';A=gallery(sA,n);
      case 303, sA='chow';A=gallery(sA,n);
      case 304, sA='circul';A=gallery(sA,n);
      case 305, sA='cycol';A=gallery(sA,n);
      case 306, sA='dramadah';A=gallery(sA,n,1);
      case 307, sA='dramadah';A=gallery(sA,n,3);
      case 308, sA='forsythe';A=gallery(sA,n);
        %   case 309, sA='invol';A=gallery(sA,n);
      case 309, sA='leslie'; A = gallery(sA,n);
      case 310, sA='leslie';A=gallery(sA,n);
      case 311, sA='normaldata';A=gallery(sA,n,10);
      case 312, sA='orthog';A=gallery(sA,n,3); %symm non herm
      case 313, sA='orthog';A=gallery(sA,n,4);
      case 314, sA='orthog';A=gallery(sA,n,-2);
      case 315, sA='randcolu';A=gallery(sA,n);
      case 316, sA='randhess';A=gallery(sA,n);
      case 317, sA='rando';A=gallery(sA,n,1);
      case 318, sA='rando';A=gallery(sA,n,2);
      case 319, sA='rando';A=gallery(sA,n,3);
      case 320, sA='randsvd';A=gallery(sA,n,1);
      case 321, sA='randsvd';A=gallery(sA,n,2);
      case 322, sA='randsvd';A=gallery(sA,n,3);
      case 323, sA='randsvd';A=gallery(sA,n,4);
      case 324, sA='randsvd';A=gallery(sA,n,5);
      case 325, sA='smoke';A=gallery(sA,n);
      case 326, sA='smoke';A=gallery(sA,n,1);
      case 327, sA='toeppen';A=full(gallery(sA,n));
      case 328, sA='uniformdata';A=gallery(sA,n,1000);

        % Singular
      case 401, sA='gearmat';A=gallery(sA,n);
      case 402, sA='neumann';A=gallery(sA,ceil(sqrt(n))^2);

    end
end