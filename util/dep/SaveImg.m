function SaveImg(requested, varargin)
%
% SaveImg(requested) saves the current matlab figure as the requested image
% formats. formats should be string(s in a cell array)
% Possible formats are 'png', 'epsc', 'fig', 'eps', 'svg', 'pdf'
% 
% Input: 
% - requested [cell array with strings] image formats to save figure with
% - savename  [string] (optional)
% 
% SaveImg({'png', 'fig', 'svg'}) saves the current figure as a high quality png,
% a matlab figure and a Scalable Vector Graphics file.
%
% SaveImg({'png', 'fig', 'epsc'}, figName) saves the current figure with 
% the given figName as the filename
%
% Leander de Kraker
% 2019-9-12
%

h = gcf;
figure(h) % bring the figure that is about to be saved to the top

if nargin == 2 % Use 4th input as savename
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
if any(strcmp(requested, 'png')) % Save the current figure as a png with a resolution of 400
    h.Renderer = 'painters'; % painters renderer does anti-aliasing
    print(h, '-dpng',[savename '.png'], '-r450')
    fprintf('saved current figure (%d) as HQ png\nin: %s\n\n', nr, savename)
end
if any(strcmp(requested, 'fig'))
    savefig(savename)
    fprintf('saved current figure (%d) as matlab figure\nin: %s\n\n', nr, savename)
end
if any(strcmp(requested, 'epsc'))
    h.Renderer = 'opengl';
    saveas(h, savename, 'epsc')
    fprintf('saved current figure (%d) as colored vector graphics epsc\nin: %s\n\n', nr, savename)
end
if any(strcmp(requested, 'eps'))
    h.Renderer = 'opengl';
    saveas(h, savename, 'eps')
    fprintf('saved current figure (%d) as grayscale vector graphics eps\nin: %s\n\n', nr, savename)
end
if any(strcmp(requested, 'svg'))
    h.Renderer = 'painters';
    saveas(h, savename, 'svg')
    fprintf('saved current figure (%d) as Scalable Vector Graphics SVG\nin: %s\n\n', nr, savename)
end
if any(strcmp(requested, 'pdf'))
    h.Renderer = 'opengl';
    saveas(h, savename, 'pdf')
    fprintf('saved current figure (%d) as PDF\nin: %s\n\n', nr, savename)
end




