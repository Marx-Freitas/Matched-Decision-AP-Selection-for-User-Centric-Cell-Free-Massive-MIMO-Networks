%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Creates access point (AP) positions randomly on any 3D plane delimited
% by a user-defined parallelogram region. The distribution of APs can vary 
% gradually from a pure random process (Binomial Point Process) to another 
% with more regularity (Hard Core Point Process), according to the user.
%
%
% SYNTAX:
%   apsPositions = aps_distribution(numberOfAps)
%   apsPositions = aps_distribution(numberOfAps, dimensions)
%   apsPositions = aps_distribution(numberOfAps, dimensions, startPosition)
%   apsPositions = aps_distribution(---, name, value)
%   apsPositions = aps_distribution(---, name1, value1, name2, value2)
%
%
% INPUT:
%   numberOfAps:
%      Define the number of APs will be generated.
%
%   dimensions:
%      Matrix of dimensions 2-by-3, in which the lines represent 
%      three-dimensional vectors. These two vectors define the plane and 
%      the parallelogram region in which the APs will be generated. Each 
%      vector represents a side of the parallelogram whose vertex lies at 
%      the origin of the coordinate system. The norms of these vectors 
%      represent the lengths of the sides of the parallelogram. For 
%      convenience, the user can enter a 2-by-2 matrix instead of 2-by-3, 
%      in which case a third null coordinate will automatically be added 
%      to the two two-dimensional vectors. The user can also enter a 
%      2-by-1 column in the form [width ; height], in this case the matrix 
%      elements represent the width and height of a rectangle in the xy 
%      plane, that is, it is equivalent to a 2-by-3 matrix in the form 
%      [width 0 0; 0 length 0]. When not specified, it assumes the default 
%      value [1 0 0 ; 0 1 0], that is, a square with side 1 in the xy
%      plane.
%
%   startPosition:
%      Vector of dimensions 1-by-3, representing the initial position of
%      the simulation parallelogram region. For convenience, the user can
%      enter a 1-by-2 vector in the form [x0 y0], which is equivalent to
%      [x0 y0 0]. The user can also enter a scalar z0, in this case it is
%      equivalent to entering [0 0 z0]. When not specified, it assumes the
%      default value [0 0 0]. 
%
%   name, value:
%      Input argument pair that specifies some settings. Can assume:
%      ---------------------------------------------------------------
%      |      name      |                   value                    |
%      ---------------------------------------------------------------
%      |  'regularity'  |      Real value in the range [0, 1].       |
%      |                                                             |
%      |   'attempts'   |        Non-negative integer value.         |
%      ---------------------------------------------------------------
%      The name 'regularity' defines the degree of regularity in the
%      distribution of the generated APs. Being 0 a distribution with
%      minimal regularity (Binomial Point Process) and 1 with maximum
%      regularity achieved by the algorithm (Hard Core Point Process). The 
%      closer to 1, the longer the processing time. When not specified, 
%      the default is 1/2. 
%      The name 'attempts' defines the maximum number of attempts used by 
%      the algorithm to satisfy the regularity condition in each AP. When
%      not specified, the default is 30. Particular attention should be
%      paid to increasing the maximum number of attempts, with the risk of
%      exponential growth in processing time. 
%
%
% OUTPUT:
%   apsPositions:
%      Matrix of dimension numberOfAps-by-3, where the rows are the
%      three-dimensional coordinates of each generated AP. 
%
%
% EXEMPLES:
% % Exemple 1 - Regularly distribute 400 APs on the roof of a 400m x 200m
% % shed and 6m ceiling height:
%      numberOfAps = 400;
%      shedDimensions = [400; 200];
%      heightOfAps = 12;
%      apsPositions = aps_distribution(numberOfAps, ...
%          shedDimensions, heightOfAps, 'regularity', 1);
% 
%      figure
%      plot3(apsPositions(:, 1), apsPositions(:, 2), ...
%          apsPositions(:, 3), '*')
%      grid
%
% % Exemple 2 - Distribute 100 APs on a 5m square panel on a side, 
% % inclined 45 degrees from the floor: 
%      numberOfAps = 100;
%      side= 5;
%      apsPositions = aps_distribution(numberOfAps, ...
%          [side 0 0; 0 -side side], [0 side 0]);
% 
%      figure
%      plot3(apsPositions(:, 1), apsPositions(:, 2), ...
%          apsPositions(:, 3), '*')
%      grid
%
% % Exemple 3 - Distribute 1000 APs on a parallelogram randomly defined in  
% % space:
%      apsPositions = aps_distribution(1000, [rand(2, 3)]);
%
%      figure
%      plot3(apsPositions(:, 1), apsPositions(:, 2), ...
%          apsPositions(:, 3), '*')
%      grid
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Created by:    Gilvan Soares Borges
% Contact:       gilvan.borges@ifpa.edu.br
% Last modified: 07/10/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function apsPositions3D = GS_apsDistribution(numberOfAps, varargin)

[width, height, linearTransformMatrix, startPosition, regularity, attempts] = ...
    check_input_arguments(numberOfAps, varargin); % Auxiliary function: Formats input arguments and avoids possible input errors.


%% DEFINING THE DEFAULT PARTITION CELL:

du = sqrt(width*height/numberOfAps);

if (width >= du) && (height >= du) % default case.
    dx = du;
    dy = du;
else % particular cases.   
    if (width  < du)
        dx = width;
        dy = height/numberOfAps;
    end
    if (height < du)
        dx = width /numberOfAps;
        dy = height;
    end
end


%% GENERATING 2D AP POSITIONS WITH A UNIFORM BPP (BINOMIAL POINT PROCESS):

xPartition = 0:dx:width;
yPartition = 0:dy:height;

nXPartition = length(xPartition);
nYPartition = length(yPartition);
initialNumberAps = nXPartition*nYPartition;

[xCellPositions, yCellPositions] = meshgrid(xPartition, yPartition);

xApsPositions = xCellPositions + dx*rand(nYPartition, nXPartition);
yApsPositions = yCellPositions + dy*rand(nYPartition, nXPartition);

apsPositions2D = [xApsPositions(:) yApsPositions(:)];


%% INCREASING FIRST-ORDER STATIONARITY BASED ON HCPP (HARD CORE POINT PROCESS):

minDist0 = regularity*du;

listNeighbors = list_neighbors(nXPartition, nYPartition); % Auxiliary function: Presents the neighbors of each AP. 

j = 0;
pairsCloseNeighbors = zeros(4*initialNumberAps - 3*(nXPartition + nYPartition) + 2, 2); % Pre-allocation with the maximum possible number of access points nearby.
for i = 1:length(listNeighbors) % Identify nearby APs.
    neighborsOfI = listNeighbors{i}(listNeighbors{i} > i);
    neighborsNearI = neighborsOfI(  dist( apsPositions2D(neighborsOfI, :), apsPositions2D(i, :).' ) < minDist0  );
    nNeighborsNearI = length(neighborsNearI);
    
    pairsCloseNeighbors(j+1:j+nNeighborsNearI, 1) = i;
    pairsCloseNeighbors(j+1:j+nNeighborsNearI, 2) = neighborsNearI;
    
    j = j + nNeighborsNearI;
end
pairsCloseNeighbors(j+1:end, :) = [];
shiftedAps = best_shifted_aps(pairsCloseNeighbors); % Auxiliary function: Among the nearby APs, those that will be repositioned are chosen.

for i = shiftedAps % Reposition the selected nearby APs.
    counter = 0;
    lessWorstDist = min(  dist( apsPositions2D(listNeighbors{i}, :), apsPositions2D(i, :).' )  );
    
    while (counter < attempts) && (lessWorstDist <= minDist0) % Reposition while you are near other APs or if the attempt limit is reached.
        newApPosition = [(xCellPositions(i) + dx*rand())  ; (yCellPositions(i) + dy*rand())];
        
        minDist = min(  dist( apsPositions2D(listNeighbors{i}, :), newApPosition )  );
        
        if (minDist > lessWorstDist)
            lessWorstDist = minDist;
            apsPositions2D(i, :) = newApPosition.';
        end
        
        counter = counter + 1; 
    end
end


%% ADJUSTING THE NUMBER OF APS:

invalidAps = [... % APs outside the rectangular region.
                   nYPartition * find(  apsPositions2D((1:nXPartition-1).*nYPartition, 2) > height  ) ; ...
              initialNumberAps - find(  apsPositions2D(initialNumberAps-(1:nYPartition-1), 1) > width  ) ; ...
              initialNumberAps + find( (apsPositions2D(initialNumberAps, 1) > width) || (apsPositions2D(initialNumberAps, 2) > height) )  -  1 ...
                                                                                                                                           ];

initialValidNumberAps = initialNumberAps - length(invalidAps);
missingAps = numberOfAps - initialValidNumberAps;

apsPositions2D = [apsPositions2D ; zeros(missingAps, 2)];
for i = 1:(missingAps) % Adding missing APs.
    lessWorstDist = 0; counter = 0;   
    
    while (lessWorstDist <= minDist0) && (counter < attempts)    
        ind = randi(initialNumberAps);    
        newApPosition = [(xCellPositions(ind) + dx*rand()) ; (yCellPositions(ind) + dy*rand())];
        
        if (  ( newApPosition(1) <= width ) && ( newApPosition(2) <= height )  )
            minDist = min(  dist( [apsPositions2D(listNeighbors{ind}, :) ; apsPositions2D(ind, :)], newApPosition )  );

            if (minDist > lessWorstDist)
                lessWorstDist = minDist;
                apsPositions2D(initialNumberAps + i, :) = newApPosition.';
            end 
        end
        
        counter = counter + 1;
    end 
end

apsPositions2D(invalidAps, :) = [];                                          % Excluding APs outside the rectangular region.
apsPositions2D(randi(initialValidNumberAps, -missingAps, 1), :) = [];        % Excluding surplus APs.


%% LINEAR TRANSFORMATION FROM 2D TO 3D:

apsPositions3D = startPosition + apsPositions2D*linearTransformMatrix;


%% FOR QUICK TEST:

% figure
% plot(apsPositions2D(:, 1), apsPositions2D(:, 2), '*')
% hold on
% plot(xCellPositions, yCellPositions, 'r')
% plot(xCellPositions.', yCellPositions.', 'r')
% axis([0, width, 0, height])
% 
% figure
% plot3(apsPositions3D(:, 1), apsPositions3D(:, 2), apsPositions3D(:, 3), '*')
% grid


end











%% ======================== AUXILIARY FUNCTIONS 1 ========================

% Formats input arguments and avoids possible input errors.
function [width, height, linearTransformMatrix, startPosition, regularity, attempts] = check_input_arguments(numberOfAps, varargin)

%% DEFAULT VALUES:

dimensions = [1 0 0; 0 1 0];  booleanDimensionsDefault = true;
startPosition = [0 0 0];      booleanP0Default         = true;
regularity = 1/2; % With "regularity = 1/2", only 15% of APs need to be repositioned. With "regularity = 1", 60% of APs need to be repositioned.
attempts   = 30;  % With up to 30 attempts, 100% of nearby APs are guaranteed to be repositioned correctly, as long as "regularity = 1/2".


%% ALLOCATING INPUT ARGUMENTS:

names = []; values = [];
numberExtraInput = length(varargin{1});

if (numberExtraInput == 1) % aps_distribution(numberOfAps, dimensions)
    dimensions = varargin{1}{1}; booleanDimensionsDefault = false;

elseif (numberExtraInput == 2) 
    if ~ischar(varargin{1}{1}) % aps_distribution(numberOfAps, dimensions, startPosition)
        dimensions       = varargin{1}{1}; booleanDimensionsDefault = false;
        startPosition = varargin{1}{2};    booleanP0Default         = false;
        
    else % aps_distribution(numberOfAps, Name, Value)
       names{1}  = varargin{1}{1};
       values{1} = varargin{1}{2};
       
    end

elseif (numberExtraInput == 3) % aps_distribution(numberOfAps, dimensions, Name, Value)
    dimensions = varargin{1}{1}; booleanDimensionsDefault = false;
    names{1}   = varargin{1}{2};
    values{1}  = varargin{1}{3};

elseif (numberExtraInput == 4)
    if ~ischar(varargin{1}{1}) % aps_distribution(numberOfAps, dimensions, startPosition, Name, Value)
        dimensions       = varargin{1}{1}; booleanDimensionsDefault = false;
        startPosition = varargin{1}{2};    booleanP0Default         = false;
        names{1}   = varargin{1}{3};
        values{1}  = varargin{1}{4};
        
    else % aps_distribution(numberOfAps, Name1, Value1, Name2, Value2)
        names{1}  = varargin{1}{1};
        values{1} = varargin{1}{2};
        names{2}  = varargin{1}{3};
        values{2} = varargin{1}{4};
    end

elseif (numberExtraInput == 5) % aps_distribution(numberOfAps, dimensions, Name1, Value1, Name2, Value2)
    dimensions = varargin{1}{1}; booleanDimensionsDefault = false;
    names{1}   = varargin{1}{2};
    values{1}  = varargin{1}{3};
    names{2}   = varargin{1}{4};
    values{2}  = varargin{1}{5};
    
elseif (numberExtraInput == 6) % aps_distribution(numberOfAps, dimensions, startPosition, Name1, Value1, Name2, Value2)
    dimensions       = varargin{1}{1}; booleanDimensionsDefault = false;
    startPosition = varargin{1}{2};    booleanP0Default         = false;
    names{1}         = varargin{1}{3};
    values{1}        = varargin{1}{4};
    names{2}         = varargin{1}{5};
    values{2}        = varargin{1}{6};
    
elseif (numberExtraInput > 6)
    error('The number of input arguments is invalid.')
end


%% CHECKING THE INPUT ARGUMENTS "numberOfAps":

if ( ~isnumeric(numberOfAps) || ~isreal(numberOfAps) || ~isscalar(numberOfAps) ...
        || (numberOfAps <= 0) || ((numberOfAps - fix(numberOfAps)) ~= 0) )
    error('The input argument "numberOfAps" is invalid.')
end


%% CHECKING THE INPUT ARGUMENT "dimensions":

if ~booleanDimensionsDefault
    if (   ~isnumeric(dimensions) || ~isreal(dimensions) || ...
            ~(size(dimensions, 1) == 2) || ~(size(dimensions, 2) < 4)   )
        error('The input argument "dimensions" is invalid.')
    end
    
    if (size(dimensions, 2) == 1)
        if ( (dimensions(1) > 0) || (dimensions(2) > 0) )
            dimensions = [dimensions(1) 0  0; 0 dimensions(2) 0];
        else
            error('The width and height dimensions of the parallelogram region of the simulation should be positive.')
        end
        
    elseif (size(dimensions, 2) == 2)
        dimensions = [dimensions [0;0]];
    end
    
    if (rank(dimensions) < 2)
        error('The sides of the simulation parallelogram region cannot be collinear.')
    end
end

width  = norm(dimensions(1, :));
height = norm(dimensions(2, :));

linearTransformMatrix = [dimensions(1, :)./width ; dimensions(2, :)./height];


%% CHECKING THE INPUT ARGUMENT "startPosition":

if ~booleanP0Default
    if (  ~isnumeric(startPosition) || ~isreal(startPosition) ...
            || ~isvector(startPosition) || (length(startPosition) > 3)  )
        error('The input argument "startPosition" is invalid.')
    end
    
    if (length(startPosition) == 1); startPosition = [0 0 startPosition]; end
    if (length(startPosition) == 2); startPosition(3) = 0;                end
end


%% CHECKING THE INPUT ARGUMENT PAIR "name" and "value":

for i = 1:length(names)
    switch names{i}
        case 'regularity'
            regularity = values{i};
            
            if (  ~isnumeric(regularity) || ~isreal(regularity) || ~isscalar(regularity)  )
                error('The input argument "regularity" is invalid.')
            end
            if (regularity < 0)
                regularity = 0;
                warning('The "regularity" input argument cannot take values ??outside the range of zero to one. Value zero was assumed.')
            end
            if (regularity > 1)
                regularity = 1;
                warning('The "regularity" input argument cannot take values ??outside the range of zero to one. Value one was assumed.')
            end
            
        case 'attempts'
            attempts = values{i};
            if ( ~isnumeric(attempts) || ~isreal(attempts) || ~isscalar(attempts) ...
                    || (attempts < 0) || ((attempts - fix(attempts)) ~= 0) )
                error('The input argument "attempts" is invalid.')
            end
            
        otherwise
            error([ 'The input argument "' names{i} '" is invalid.' ])
    end
end


end





%% ======================== AUXILIARY FUNCTIONS 2 ========================

% For each AP, it determines its neighbors.
function listNeighbors = list_neighbors(nXPartition, nYPartition)

numberAps = nXPartition*nYPartition;
listNeighbors = cell(numberAps, 1);

for k = 1:numberAps
    if (k == 1)
        listNeighbors{k} = [(k+1) ; ...
                            (k+nYPartition) ; ...
                            (k+nYPartition+1)    ];
        
    elseif(k == nYPartition)
        listNeighbors{k} = [(k-1); ...
                            (k+nYPartition-1); ...
                            (k+nYPartition)       ];
        
    elseif (k == numberAps-nYPartition+1)
        listNeighbors{k} = [(k-nYPartition); ...
                            (k-nYPartition+1); ...
                            (k+1)                    ];
        
    elseif (k == numberAps)
        listNeighbors{k} = [(k-nYPartition-1); ...
                            (k-nYPartition); ...
                            (k-1)                    ];
        
    elseif ( (k > 1) && (k < nYPartition) )
        listNeighbors{k} = [(k-1); ...
                            (k+1); ...
                            (k+nYPartition-1); ...
                            (k+nYPartition); ...
                            (k+nYPartition+1)     ];
        
    elseif (mod(k-1, nYPartition) == 0)
        listNeighbors{k} = [(k-nYPartition); ...
                            (k-nYPartition+1);  ...
                            (k+1); ...
                            (k+nYPartition); ...
                            (k+nYPartition+1)     ];

    elseif (mod(k, nYPartition) == 0)
        listNeighbors{k} = [(k-nYPartition-1); ...
                            (k-nYPartition); ...
                            (k-1); ...
                            (k+nYPartition-1); ...
                            (k+nYPartition)       ];
        
    elseif ( (k > numberAps-nYPartition+1) && (k < numberAps) )
        listNeighbors{k} = [(k-nYPartition-1); ...
                            (k-nYPartition); ...
                            (k-nYPartition+1); ...
                            (k-1); ...
                            (k+1)                 ];
        
    else
        listNeighbors{k} = [(k-nYPartition-1); ...
                            (k-nYPartition); ...
                            (k-nYPartition+1); ...
                            (k-1); ...
                            (k+1); ...
                            (k+nYPartition-1); ...
                            (k+nYPartition); ...
                            (k+nYPartition+1)     ];
    end
    
end


end





%% ======================== AUXILIARY FUNCTIONS 3 ========================

% Among the nearby access points, it is chosen which are the most suitable 
% for relocating positions. The basic idea is to prioritize the indexes 
% that are most repeated.
function shiftedAps = best_shifted_aps(pairsCloseNeighbors)

if isempty(pairsCloseNeighbors)
    shiftedAps = [];
    return
end

frequency =  histcounts(pairsCloseNeighbors(:), (0:max(max(pairsCloseNeighbors))) + 1/2).';
pairsCloseNeighbors = sortrows([ max(frequency(pairsCloseNeighbors(:, 1)), frequency(pairsCloseNeighbors(:, 2))) pairsCloseNeighbors ], 1, 'descend');
pairsCloseNeighbors(:, 1) = [];

nPairsCloseNeighbors = size(pairsCloseNeighbors, 1);
shiftedAps = zeros(1, nPairsCloseNeighbors);
for i = 1:nPairsCloseNeighbors
    Ap1Priority = frequency(pairsCloseNeighbors(i, 1));
    Ap2Priority = frequency(pairsCloseNeighbors(i, 2));
    
    if (Ap1Priority == Ap2Priority)
        if (Ap1Priority == 1)
            shiftedAps(i) = pairsCloseNeighbors(i, randi(2));
        else
            balance1 = length(find(shiftedAps(1:i-1) == pairsCloseNeighbors(i, 1))) - length(find(pairsCloseNeighbors(1:i-1, :) == pairsCloseNeighbors(i, 1)));
            balance2 = length(find(shiftedAps(1:i-1) == pairsCloseNeighbors(i, 2))) - length(find(pairsCloseNeighbors(1:i-1, :) == pairsCloseNeighbors(i, 2)));
            
            if (balance1 > balance2)
                shiftedAps(i) = pairsCloseNeighbors(i, 1);
            else
                shiftedAps(i) = pairsCloseNeighbors(i, 2);
            end
            
        end
        
    elseif (Ap1Priority > Ap2Priority)
        shiftedAps(i) = pairsCloseNeighbors(i, 1);
    else
        shiftedAps(i) = pairsCloseNeighbors(i, 2);
    end
    
end

shiftedAps = unique(shiftedAps);

end