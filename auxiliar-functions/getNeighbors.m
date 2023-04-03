function n = getNeighbors(gamma, w, h, olapx, olapy, dx, dy)
% [Description]
ovlapx = olapx*w/100; ovlapy = olapy*h/100;
%
%if nargin < 4
n = cell(8,1);
n{1} = gamma' + (w-ovlapx)*dx +  (h-ovlapy)*dy; % northeast
n{2} = gamma' + (w-ovlapx)*dx ;                 % east
n{3} = gamma' + (w-ovlapx)*dx + (-h+ovlapy)*dy; % southeast
n{4} = gamma' + (h-ovlapy)*dy;                  % north
n{5} = gamma' + (-h+ovlapy)*dy;                 % south
n{6} = gamma' + (-w+ovlapx)*dx + (h-ovlapy)*dy; % northwest
n{7} = gamma' + (-w+ovlapx)*dx;                 % west
n{8} = gamma' + (-w+ovlapx)*dx +(-h+ovlapy)*dy; % southwest
%else
%     switch varargin{4}
%         case 'cardinal'
%             aux{1} = map{indrow - 1, indcol    }; % north
%             aux{2} = map{indrow    , indcol + 1}; % east
%             aux{3} = map{indrow + 1, indcol    }; % south
%             aux{4} = map{indrow    , indcol - 1}; % west
%         case 'diagonal'
%             aux{1} = map{indrow - 1, indcol - 1}; % northwest
%             aux{2} = map{indrow - 1, indcol + 1}; % northeast
%             aux{3} = map{indrow + 1, indcol + 1}; % southeast
%             aux{4} = map{indrow + 1, indcol - 1}; % southwest
%     end
% end

%%
% n = {};
% indrow = varargin{1};
% indcol = varargin{2};
% aux_map = varargin{3};
% map = cell(size(aux_map,1) + 2, size(aux_map,2) + 2);
% for i=1:size(map,1)
%     for j=1:size(map,2)
%         map{i,j} = [NaN NaN];
%     end
% end
% map(2:end-1,2:end-1) = aux_map;
% indrow = indrow + 1;
% indcol = indcol + 1;
% if nargin < 4
%     n = cell(8,1);
%     if ~isnan(map{indrow, indcol})
%         n{1} = map{indrow - 1, indcol + 1}; % northeast
%         n{2} = map{indrow    , indcol + 1}; % east
%         n{3} = map{indrow + 1, indcol + 1}; % southeast
%         n{4} = map{indrow - 1, indcol    }; % north
%         n{5} = map{indrow + 1, indcol    }; % south
%         n{6} = map{indrow - 1, indcol - 1}; % northwest
%         n{7} = map{indrow    , indcol - 1}; % west
%         n{8} = map{indrow + 1, indcol - 1}; % southwest
%     end
% else
%     switch varargin{4}
%         case 'cardinal'
%             n{1} = map{indrow - 1, indcol    }; % north
%             n{2} = map{indrow    , indcol + 1}; % east
%             n{3} = map{indrow + 1, indcol    }; % south
%             n{4} = map{indrow    , indcol - 1}; % west
%         case 'diagonal'
%             n{1} = map{indrow - 1, indcol - 1}; % northwest
%             n{2} = map{indrow - 1, indcol + 1}; % northeast
%             n{3} = map{indrow + 1, indcol + 1}; % southeast
%             n{4} = map{indrow + 1, indcol - 1}; % southwest
%     end
% end
% 
% for i=1:length(n)
%     if ~isnan(n{i})
%         n{end + 1} = n{i};
%     end
% end

end