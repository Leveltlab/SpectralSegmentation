function colmap = cmapL(colors, n)
% 
% Create a custom colormap from multiple RGB triplets or a string
% If only one RGB triplet is given, that color will be the maximum color, 
% the minimum color will be black
% 
% colmap = cmapL(colors, n);
%
% Input:
%   colors [y x 3] 2D double: the RGB triples that are the colormap
%   n (scalar, or [1 x (y-1)]): how many colors to output. 
%           In the case of vector input: how many colors for each given
%           color transition.
%
%
% EXAMPLES:
% 
% img = zeros(124); % first create an image to use the colormap on
% img([63:65], [63:65]) = 1;
% img = imgaussfilt(img,33);
% 
% n = 12; % 12 colors for the colormap
% colors = [0 1 0]; % highest values will be green. standard low = black
% colmap = cmapL(colors, n);
% figure; surf(img,'edgecolor','none'); colormap(colmap); colorbar
% title("cmapL([0 1 0],12)")
%
% n = 256; % 256 colors for the colormap
% colors = [1 0 0; 1 1 0; 0 1 1; 0 0 1] % creating jetlike colormap
% colmap = cmapL(colors, n);
% figure; surf(img,'edgecolor','none'); colormap(colmap); colorbar
% title("cmapL( [1 0 0; 1 1 0; 0 1 1; 0 0 1] ,256)")
%
% n = 256; % 256 colors for the colormap
% colors = 'b'; % ask for blue color map
% colmap = cmapL(colors, n);
% figure; surf(img,'edgecolor','none'); colormap(colmap); colorbar
% title("cmapL('b',256)")
%
% colmap = cmapL('dutch', 3); % ask for dutch flag color map, only 3 colors
% figure; surf(img,'edgecolor','none'); colormap(colmap); colorbar
% title("cmapL('dutch',3)")
%
% n = [20, 236]; % Using 20 colors for red-white, 236 colors for white-blue
% colmap = cmapL('dutch', n);
% figure; surf(img,'edgecolor','none'); colormap(colmap); colorbar
% title("cmapL('dutch', [20, 236]): 20 for red-white, 236 for white-blue")
%
%
%
% possible colors by string request:
% w  = white
% r  = red
% g  = green
% b  = blue, but extra light
% y  = yellow
% c  = cyan (blue and green)
% m  = magenta
% dutch     = official reddish, white, blueish
% italian   = official reddish, white, greenish
% italian roast = from yellowwhite, red  , black, green, cyanwhite
% blue roast = from red, black, blue(white) safer for colorblind
% viridis   = purplish to green to yellow, from matplotlib
% pastel    = pastelish reddish, to blueish
% inferno   = black to purple to yellow, standard matplotlib
% greenFancy = from black to blue to green to white
% greenFancyDark = darker version of greenFancy
% painbow   = worst color scale award: https://xkcd.com/2537/
% moonbow   = colors for Pink Floyd's Dark Side Of The Moon rainbow
% safe      = suitable for colorblind in 11 steps
%
% Leander de Kraker
% 2018-10-22, edited 2020-5-1
% 2022-2-11: added ability to use different amounts of colors per color
%

% Change string request to
if ischar(colors)
    switch lower(colors)
        case 'w'
            colors = [1 1 1; 0 0 0];
        case 'r'
            colors = [1 0.1 0.1; 0 0 0];
        case 'g'
            colors = [0 1 0; 0 0 0];
        case 'b'
            colors = [0.4 0.4 1; 0 0 0]; 
        case 'y'
            colors = [1 1 0; 0 0 0];
        case 'c'
            colors = [0 1 1; 0 0 0];
        case 'm'
            colors = [1 0 1; 0 0 0];
        case 'dutch' % dutch flag colors
            colors = [0.68 0.1 0.15; 1 1 1; 0.12 0.27 0.545];
        case 'italian' % italian flag color
            colors = [0 0.57 0.27; 1 1 1; 0.8 0.16 0.17];
        case 'italian roast' % italian flag, but more extreme and black instead of white
            colors = [0.8 1 1; 0.6 1 0.6; 0 0.8 0; 0 0 0; 0.8 0 0; 1 0.6 0.6; 1 1 0.8]; 
        case 'blue roast' % inspired by italian roast, but safe for colorblind
            colors = [1 1 1; 0.7 0.8 1; 0.3 0.3 1; 0 0 0; 0.8 0 0; 1 0.6 0.6; 1 1 0.8];
        case 'safe' % Colorblind proof, line suitable
            colors = [51,  34,  136;...
                      102, 153, 204;...
                      136, 204, 238;...
                      68,  170, 153;...
                      17,  119, 51;...
                      153, 153, 51;...
                      221, 204, 119;...
                      204, 102, 119;...
                      170, 68,  153;...
                      136, 34,  85;...
                      102, 17,  0]./255;
        case 'pastel'
            colors = [1 0.6 0.56; 0.56 1 1];
        case 'viridis' % A slightly darkened version of viridis
            colors = [0.9961    0.9094    0.15;...
                      0.7087    0.8740    0.18;...
                      0.2087    0.7205    0.48;...
                      0.1496    0.5118    0.56;...
                      0.22      0.25      0.4;...
                      0.1       0         0.15];
        case 'inferno'
            colors = [249 244 169;...
                      250 205 55;...
                      247 143 30;...
                      228 91 48;...
                      186 56 85;...
                      138 35 106;...
                      85 37 105;...
                      30 19 73;...
                      0 0 0]./254;
        case 'greenfancy'
            colors = [1     1    1;...
                      1     1    0.75;...
                      0.75  1    0.5;...
                      0.25  0.75 0.5;...
                      0     0.25 0.5;...
                      0     0    0];
        case 'greenfancydark'
            colors = [1     1    1;
                      0.75  1    0.5;...
                      0     0.25 0.5;...
                      0     0    0];
        case 'painbow'
            colors = [255 255 255;...
                      193 241 104;...
                      64  217 50;...
                      12  178 22;...
                      62  134 10;...
                      193 79  1;...
                      230 10  5;...
                      255 0   0;...
                      13  125 80;...
                      34  46 100;...
                      0   0   0;...
                      255 255 255]./255;
        case 'moonbow'
            colors = [121, 102, 150;...
                      86,  158, 193;...
                      142, 165, 63;...
                      231, 207, 36;...
                      209, 128, 45;...
                      194, 71,  72]./255;
            % black: [20, 23, 30]. white: [233, 226, 210]. gray [131, 144, 152]
        otherwise
            warning('unknown colormap inputted: defaulting to white')
            colors = [1 1 1; 0 0 0];
    end
end

colmap = zeros(sum(n),3); % preallocate colormap
if size(colors, 1) == 1 % Add the standard black color if necessary
    colors(2,:) = [0 0 0];
end
ncolors = size(colors, 1); % The number of different RGB colors in the colormap

lims = floor(linspace(1, sum(n), ncolors));

% Use a certain amount of colors between each two requested colors
if length(n)>1
    if length(n) == (ncolors-1) % Create the requested limits
        lims = cumsum([1 n]);
    elseif length(n) > (ncolors-1)
        warning(sprintf('Remove %d number of colors to output numbers. Or use only 1 value',...
                        length(n)-(ncolors-1)))
    elseif (length(n)-1) < ncolors
        warning(sprintf('Add %d number of colors to output numbers. Or use only 1 value',...
                        (ncolors-1)-length(n)))
    end
end

% Fill the colormap
for i = (ncolors-1):-1:1
    for c = 1:3
        % Use linear spacing between the colorlimits for each color
        colmap(lims(i):lims(i+1),c) = linspace(colors(i,c), colors(i+1,c), lims(i+1)-lims(i)+1);
    end
end

colmap = flipud(colmap);

end