function [corner, check] = GetRoiCorners(Mask, PP, decreSize, checkSize, varargin)
% Get coordinates of corners of ROIs.
%
% which ROIs? Well.. This can be selected based on size, or on ROI number
% If you want to select on size, call the function with 'size' as 5th input
% argument, and the 4th input will be minimum size to be analyzed
% If there is no 5th input argument the 4th input argument should be the
% ROI number(s) of the ROIs to check, in a vector
%
%
% I'll pick the most extreme coordinates in a ROI to start checking 
% correlation with the rest of the ROI. But the points shouldn't be 
% outside of neurons. To  ensure this, I decrease the size of the masks,
% using my BufferMask function in a reverse way


rMask = Mask; % reverse valued Mask (0 -> 1), (>=1 -> 0)
rMask(rMask>0) = 1;
rMask = 1-rMask;
[~,~,rMask] = BufferMask(rMask,decreSize); % Finding inner pixels on ROI boundries
uMask = Mask; % useMask
uMask(rMask>0)=0;

if nargin > 4
    if strcmpi(varargin{1}, 'size')
        % When ROI is bigger than given pixels in area, check it for splitting
        check = find(PP.A > checkSize(1));
        if length(checkSize) > 1
            warning('when "size" selection modus is enabled, input only minimum size\n')
        end
    else % for some strange reason the 5th input argument was not "size"
        % Check the requested neuron(s)
        check = checkSize;
    end
else
    % Check the requested neuron(s)
    check = checkSize;
end

corner = zeros(2,4,length(check));

% Get coordinates to start the correlation from
for i = 1:length(check)
    [y, x] = find(uMask==check(i));
    [corner(1,1,i), idx] = min(x);
    corner(2,1,i)        = y(idx);
    [corner(1,2,i), idx] = max(x);
    corner(2,2,i)        = y(idx);
    [corner(2,3,i), idx] = min(y);
    corner(1,3,i)        = x(idx);
    [corner(2,4,i), idx] = max(y);
    corner(1,4,i)        = x(idx);
end

end

