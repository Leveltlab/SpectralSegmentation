function score = OverlapScore(inRoi, linkMat, i)
    % Calculate the average % overlap between ROIs, 
    %
    % input:
    %       inRoi: cell array with cell arrays with overlap information
    %       linkMat: Matrix with ROI numbers, for row i a score is calculated
    %       i: digid. which row to calculate overlap score for
    %
    % iutput:
    %       score: % overlap between the ROIs in row i of the linkMat
    %
    % Leander de Kraker
    %

    score = 0;
    counter = 0;
    present = find(linkMat(i,:));
    if length(present) > 1 % If there are matches to be searched
        for m = present(1:end-1)
            for c = present(present>m)
                % Check overlap from roi i in recording m -> linked roi in recording c
                found = inRoi{m}{linkMat(i,m),c};
                idx = found(:,1)==linkMat(i,c);
                if ~isempty(idx) && any(idx) % Maybe the ROI wasn't in the ROI, abort in that case
                    try
                        score = score +  found(idx, 3);
                    catch
                        fprintf('Doubt for match %d: no overlap between %d and %d\n', i, m, c)
                    end
                else
                    fprintf('Doubt for match %d: no overlap for neuron %d and %d\n', i, m, c)
                end

                % Check overlap from roi i in recording c -> linked roi in
                % recording m (which is the roi i in recording m)
                found = inRoi{c}{linkMat(i,c),m};
                idx = found(:,1)==linkMat(i,m);
                if ~isempty(idx) && any(idx) % Maybe the ROI wasn't in the ROI, abort in that case
                    try
                        score = score +  found(idx, 3);
                    catch
                        fprintf('Doubt for match %d: no overlap between %d and %d\n', i, c, m)
                    end
                end

                counter = counter + 2;

            end
        end
        score = score./counter; % calculate average overlap
    end
        
end