function mp_semilogy(x, y, axis_lim, n_yticks, varargin)
% MP_SEMILOGY Handle values smaller than 1e-300 with SEMILOGY.

% Copyright (c) 2016-2018, Massimiliano Fasi
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
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


    plot(x, log10(y), varargin{:})
    if (isempty(axis_lim))
        axis([min(x)-1, max(x)+1,...
              double(min(log10(y(y>0))))-1, double(max(log10(y(y>0))))])
    else
        axis_lim(3:4) = log10(axis_lim(3:4));
        axis_lim = double(axis_lim);
        axis(axis_lim);
        if (n_yticks > 0)
            step = (axis_lim(4) - axis_lim(3)) / (n_yticks - 1);
            set(gca,'Ytick', axis_lim(3):step:axis_lim(4))
        end
    end
    ax = gca;
    y = get(ax, 'YTick');
    strings = cell(length(y), 1);
    for i = 1 : length(y)
        if(round(y(i)) == y(i))
            strings{i} = sprintf('10^{%d}', y(i));
        else
            strings{i} = sprintf('10^{%.1f}', y(i));
        end
    end
    set(ax,'YTickLabel',strings)
end
