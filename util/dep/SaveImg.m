function SaveImg(requested, varargin)
%
% SaveImg(requested) saves the current matlab figure as the requested image
% formats. formats should be string(s in a cell array)
% Possible formats are 'png', 'epsc', 'fig', 'eps', 'svg', 'pdf'
% 
% Input: 
% - requested [cell array with strings] image formats to save figure with
% - savename  [string] (optional)
% - resolution [scalar value] (optional) Only applies for png images
% 
% SaveImg({'png', 'fig', 'svg'}) saves the current figure as a high quality png,
% a matlab figure and a Scalable Vector Graphics file.
%
% SaveImg({'png', 'fig', 'epsc'}, figName, 500) saves the current figure with 
% the given figName as the filename, with a good resolution of 500 for png
%
% Leander de Kraker
% 2019-9-12
% 2022-8-17: added resolution as optional input.
%

h = gcf;
figure(h) % bring the figure that is about to be saved to the top

if nargin == 3
    resolution = ['-r' num2str(varargin{2})];
else
    resolution = '-r450';
end
if nargin >= 2 % Use 4th input as savename
    savename = varargin{1};
elseif nargin == 1 % Ask for name that the picture should have
    [savefilename, savepathname] = uiputfile('*');
    if savefilename==0
        fprintf('No file selected, not saving picture\n')
        return
    end
    savefilename = strsplit(savefilename, '.');
    savefilename = savefilename{1};
    savename = [savepathname savefilename];
else
    fprintf('Incorrect input, please read:\n')
    help SaveImg
    return
end

nr = get(gcf, 'Number');
fprintf('\n')
% Save the current figure
if any(strcmp(requested, 'png'))
    h.Renderer = 'painters'; % painters renderer does anti-aliasing
    print(h, '-dpng',[savename '.png'], resolution)
    fprintf('saved current figure (%d) as HQ png\nin: %s\n\n', nr, savename)
end
if any(strcmp(requested, 'jpg'))
    saveas(h, savename, 'jpg')
    fprintf('saved current figure (%d) as shit jpg\nin: %s\n\n', nr, savename)
end
if any(strcmp(requested, 'fig'))
    savefig(savename)
    fprintf('saved current figure (%d) as matlab figure\nin: %s\n\n', nr, savename)
end
if any(strcmp(requested, 'epsc'))
    h.Renderer = 'opengl';
    saveas(h, [savename, 'epsc'], 'epsc')
    fprintf('saved current figure (%d) as colored vector graphics epsc\nin: %s\n\n', nr, savename)
end
if any(strcmp(requested, 'eps'))
    h.Renderer = 'opengl';
    saveas(h, [savename, 'eps'], 'eps')
    fprintf('saved current figure (%d) as grayscale vector graphics eps\nin: %s\n\n', nr, savename)
end
if any(strcmp(requested, 'svg'))
    h.Renderer = 'painters';
    saveas(h, [savename, '.svg'], 'svg')
    fprintf('saved current figure (%d) as Scalable Vector Graphics SVG\nin: %s\n\n', nr, savename)
end
if any(strcmp(requested, 'pdf'))
    h.Renderer = 'opengl';
    saveas(h, [savename, 'pdf'], 'pdf')
    fprintf('saved current figure (%d) as PDF\nin: %s\n\n', nr, savename)
end



