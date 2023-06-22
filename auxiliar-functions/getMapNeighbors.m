function n = getMapNeighbors(varargin)
% [Description]

%%
n = {};
indrow = varargin{1};
indcol = varargin{2};
map = varargin{3};
% aux_map = varargin{3};
% map = cell(size(aux_map,1) + 2, size(aux_map,2) + 2);
% for i=1:size(map,1)
%     for j=1:size(map,2)
%         map{i,j} = [NaN NaN];
%     end
% end
% map(2:end-1,2:end-1) = aux_map;
if nargin < 4
    aux_n = cell(8,1);
    if ~any(isnan(map{indrow, indcol})) && ~isempty(map{indrow, indcol})
        aux_n{1} = map{indrow - 1, indcol + 1}; % northeast
        aux_n{2} = map{indrow    , indcol + 1}; % east
        aux_n{3} = map{indrow + 1, indcol + 1}; % southeast
        aux_n{4} = map{indrow - 1, indcol    }; % north
        aux_n{5} = map{indrow + 1, indcol    }; % south
        aux_n{6} = map{indrow - 1, indcol - 1}; % northwest
        aux_n{7} = map{indrow    , indcol - 1}; % west
        aux_n{8} = map{indrow + 1, indcol - 1}; % southwest
    end
else
    switch varargin{4}
        case 'cardinal'
            aux_n{1} = map{indrow - 1, indcol    }; % north
            aux_n{2} = map{indrow    , indcol + 1}; % east
            aux_n{3} = map{indrow + 1, indcol    }; % south
            aux_n{4} = map{indrow    , indcol - 1}; % west
        case 'diagonal'
            aux_n{1} = map{indrow - 1, indcol - 1}; % northwest
            aux_n{2} = map{indrow - 1, indcol + 1}; % northeast
            aux_n{3} = map{indrow + 1, indcol + 1}; % southeast
            aux_n{4} = map{indrow + 1, indcol - 1}; % southwest
    end
end

for i=1:length(aux_n)
    if ~any(isnan(aux_n{i})) && ~isempty(aux_n{i})
        n{end + 1} = aux_n{i};
    end
end

end