function n = getNeighbors(varargin)
% [Description]

%%
n = {};
indrow = varargin{1};
indcol = varargin{2};
aux_map = varargin{3};
map = cell(size(aux_map,1) + 2, size(aux_map,2) + 2);
for i=1:size(map,1)
    for j=1:size(map,2)
        map{i,j} = [NaN NaN];
    end
end
map(2:end-1,2:end-1) = aux_map;
indrow = indrow + 1;
indcol = indcol + 1;
if nargin < 4
    aux = cell(8,1);
    if ~isnan(map{indrow, indcol})
        aux{1} = map{indrow - 1, indcol + 1}; % northeast
        aux{2} = map{indrow    , indcol + 1}; % east
        aux{3} = map{indrow + 1, indcol + 1}; % southeast
        aux{4} = map{indrow - 1, indcol    }; % north
        aux{5} = map{indrow + 1, indcol    }; % south
        aux{6} = map{indrow - 1, indcol - 1}; % northwest
        aux{7} = map{indrow    , indcol - 1}; % west
        aux{8} = map{indrow + 1, indcol - 1}; % southwest
    end
else
    switch varargin{4}
        case 'cardinal'
            aux{1} = map{indrow - 1, indcol    }; % north
            aux{2} = map{indrow    , indcol + 1}; % east
            aux{3} = map{indrow + 1, indcol    }; % south
            aux{4} = map{indrow    , indcol - 1}; % west
        case 'diagonal'
            aux{1} = map{indrow - 1, indcol - 1}; % northwest
            aux{2} = map{indrow - 1, indcol + 1}; % northeast
            aux{3} = map{indrow + 1, indcol + 1}; % southeast
            aux{4} = map{indrow + 1, indcol - 1}; % southwest
    end
end

for i=1:length(aux)
    if ~isnan(aux{i})
        n{end + 1} = aux{i};
    end
end

end