function w = triang_win (n)
% 
% function w = triang_win (n) 
% It calculates the n points triangular window.
% 
% n   :   window length.
% 
% w   :   triangular window.

if rem(n,2)
    % It's an odd length sequence
    w = 2*(1:(n+1)/2)/(n+1);
    w = [w w((n-1)/2:-1:1)]';
else
    % It's even
    w = (2*(1:(n+1)/2)-1)/n;
    w = [w w(n/2:-1:1)]';
end