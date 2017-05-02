function [ loc ] = ApplyVarargin(loc, varargin)
%
% ApplyVarargin (loc, varargin)
%
% loc - a local structure which varargin will be coppied into
% varargin  - has form {(opt) structure, 'fieldname', value, ....}           
%                                                                        
% The input structure must has all the fieldnames supplied in varargin    
%

% read the acceptable local names
optionNames = fieldnames(loc);
varargin=varargin{1};

% arguments must come loc arragenment, ((opt) structure, "then 0+ pairs of " name', value ... )
nArgs = length(varargin);

if round(nArgs/2)~=nArgs/2
    inStruct = varargin{1};
    % Is the lenght is odd, check if the first arguement is a structure
    if(isstruct(inStruct))
        varargin = varargin(2:end);
        inNames = fieldnames(inStruct);
        
        % Copy loc values from external structure
        for i = 1:length(inNames)
            if any(strcmp(inNames{i},optionNames))
                loc.(inNames{i}) = inStruct.(inNames{i});
            end
        end
        
    else
        error('EXAMPLE needs ((opt) struct, (opt) ''propertyName'', propertyValue, ...) loc pairs).')
    end
end


% Populate from input and make sure consistent - overwite anything in 'in'
for pair = reshape(varargin,2,[]) %# pair is {name;value}
    inpName = pair{1};
    
    if any(strcmp(inpName,optionNames))
        loc.(inpName) = pair{2};
    else
        error('%s is not recognized parameter name',inpName)
    end
end


end

