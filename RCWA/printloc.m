function success = printloc(varargin)
success = true;
try
if (nargin == 1)
    loc = varargin{1};
    filename = 'loc.dat';
elseif (nargin == 2)
    loc = varargin{1};
    filename = varargin{2};
else
    loc = DefaultLoc;
    filename = 'loc.dat';
end

fn = fieldnames(loc);

fileID = fopen(filename,'wt');

for i = 1:length(fn)
fprintf(fileID, sprintf('%s: %s\n',fn{i},num2str(loc.(fn{i}))));
end

fclose(fileID);

catch
    success = false;
    try
        fclose(fileID);
    catch
        sprintf('Error and no file to close');
    end
    
    sprintf('Error, but file closed');
end

end

