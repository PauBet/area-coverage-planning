function grid = checkTaboo3(grid, ind_row, ind_col, dir1, dir2)

% Check if there are any taboo tiles that should be deleted... (at the
% expense of potential uncovered area)
% This matrix points at the taboo tiles of grid (taboo_idx(i, j) = 1 if
% there is a taboo tile at element [i, j] of the grid)
taboo_idx = zeros(size(grid));

switch dir1
    case {'north', 'south'} % horizontal sweep

        if strcmp(dir1, 'south') % spacecraft is towards roi's
            % bottom
            if ind_row > 1 && any(cellfun(@any, ...
                    grid(1:(ind_row - 1), :)), 'all')

                % Mark the taboo tile at the taboo_idx matrix
                [row, col] = find(~cellfun(@isempty, ...
                    grid(1:(ind_row - 1), :)));
                taboo_idx(row, col) = 1;
            end
        else % spacecraft is towards roi's top
            if ind_row < size(grid, 1) && any(cellfun(@any, ...
                    grid((ind_row + 1):end, :)), 'all')

                % Mark the taboo tile at the taboo_idx matrix
                [row, col] = find(~cellfun(@isempty, ...
                    grid((ind_row + 1):end, :)));
                row = row + ind_row; % adjust row indices to be
                % relative to the entire grid
                taboo_idx(row, col) = 1;
            end
        end

        if isequal(dir2, 'east') % tour is moving to the right (left -> right dir.)
            if ind_col > 1
                if isequal(dir1, 'south') && any(cellfun(@any, ...
                        grid(1:ind_row, 1:(ind_col - 1))), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid(1:ind_row, 1:(ind_col - 1))));
                    taboo_idx(row, col) = 1;
                elseif isequal(dir1, 'north') && any(cellfun(@any, ...
                        grid(ind_row:end, 1:(ind_col - 1))), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid(ind_row:end, 1:(ind_col - 1))));
                    row = row + ind_row - 1; % adjust row indices to be
                    % relative to the entire grid
                    taboo_idx(row, col) = 1;
                end
            end
        else % tour is moving to the left (right -> left dir.)
            if ind_col < size(grid, 2)
                if isequal(dir1, 'south') && any(cellfun(@any, ...
                        grid(1:ind_row, (ind_col + 1):end)), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid(1:ind_row, (ind_col + 1):end)));
                    col = col + ind_col; % adjust column indices to
                    % be relative to the entire grid
                    taboo_idx(row, col) = 1;
                elseif isequal(dir1, 'north') && any(cellfun(@any, ...
                        grid(ind_row:end, (ind_col + 1):end)), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid(ind_row:end, (ind_col + 1):end)));
                    col = col + ind_col; % adjust column indices to
                    % be relative to the entire grid
                    row = row + ind_row - 1; % adjust row indices
                    % to be relative to the entire grid
                    taboo_idx(row, col) = 1;
                end

            end
        end

    case {'east', 'west'} % vertical sweep

        if strcmp(dir1, 'east')
            if ind_col > 1 && any(cellfun(@any, ...
                    grid(:, 1:(ind_col - 1))), 'all')

                % Mark the taboo tile at the taboo_idx matrix
                [row, col] = find(~cellfun(@isempty, ...
                    grid(:, 1:(ind_col - 1))));
                taboo_idx(row, col) = 1;
            end
        else
            if ind_col < size(grid, 2) && any(cellfun(@any, ...
                    grid(:, (ind_col + 1):end)), 'all')

                % Mark the taboo tile at the taboo_idx matrix
                [row, col] = find(~cellfun(@isempty, ...
                    grid(:, (ind_col + 1):end)));
                col = col + ind_col; % adjust column indices to
                % be relative to the entire grid
                taboo_idx(row, col) = 1;
            end
        end

        if isequal(dir2, 'south') % downsweep = true. tour is moving to the bottom
            if ind_row > 1
                if isequal(dir1, 'east') && ...
                        any(cellfun(@any, grid(1:(ind_row - 1), ...
                        1:ind_col)), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid(1:(ind_row - 1), 1:ind_col)));
                    taboo_idx(row, col) = 1;

                elseif isequal(dir1, 'west') && ...
                        any(cellfun(@any, grid(1:(ind_row - 1), ...
                        ind_col:end)), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid(1:(ind_row - 1), ind_col:end)));
                    col = col + ind_col - 1;
                    taboo_idx(row, col) = 1;
                end
            end
        else % tour is moving to the top
            if ind_row < size(grid, 1)

                if isequal(dir1, 'east') && ...
                        any(cellfun(@any, grid((ind_row + 1):end, ...
                        1:ind_col)), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid((ind_row + 1):end, 1:ind_col)));
                    row = ind_row + row;
                    taboo_idx(row, col) = 1;
                elseif isequal(dir1, 'west') && ...
                        any(cellfun(@any, grid((ind_row + 1):end, ...
                        ind_col:end)), 'all')

                    % Mark the taboo tile at the taboo_idx matrix
                    [row, col] = find(~cellfun(@isempty, ...
                        grid((ind_row + 1):end, ind_col:end)));
                    row = ind_row + row;
                    col = ind_col + col - 1;
                    taboo_idx(row, col) = 1;
                end
            end
        end
end

% Delete taboo tiles
grid(taboo_idx == 1) = {[]};

% After deleting the taboo tiles, check if there are empty rows and/or
% columns and erase those from the grid...
empty_rows    = all(cellfun('isempty', grid), 2);
empty_columns = all(cellfun('isempty', grid), 1);
grid(empty_rows, :) = []; grid(:, empty_columns) = [];
end