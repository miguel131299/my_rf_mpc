function alpha_int = bz_int(alpha,x0,s_max)
% performs numerical integration using the trapezoidal rule.

% INPUT.
% alpha: input signal
% x0: initial condition
% s_max: integration range

%% Input validation
if nargin == 2
    s_max = 1;
end

[n,m] = size(alpha);        % make sure alpha is a row vector
if n > m
    alpha = alpha';         % transpose if necessary
end

M = length(alpha);
AA = zeros(M+1,M+1);

%% Trapezoidal Rule 
for ii = 1:M
    AA(ii,ii:ii+1) = [-1 1];       % fill with trapezoidal rule coefficients
end

AA = M/s_max*AA;    % scale coefficient matrix
AA(M+1,1) = 1;      % set boundary condition (integrated signal starts with
                    % initial condition x0

% numerical integration
% AA * [alpha X0]' = alpha int
% alpha_int - integrated signal
alpha_int = (AA\[alpha x0]')';  % linear system solved bz matrix inversion


