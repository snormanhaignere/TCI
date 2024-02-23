function [I, C, C_value, all_keys, paramstring] = parse_optInputs_keyvalue(optargs, I, varargin)

% Parse optional inputs specified as a key-value pair. Key must be a string.
% 
% I is an optional input argument with a structure containing the default values
% of the parameters, which are overwritten by the optional arguments. If I is
% specified then every key must match one of the fields of the structure I. 
% 
% C is a structure that tells you whether each field has been modified.
% 
% C_value has the value for all of the fields that have been changed
% 
% -- Example --
% 
% % parse key value pairs
% optargs = {'key1', {1 2 3}, 'key2', 1:3, 'key3', 'abc'};
% [I, C] = parse_optInputs_keyvalue(optargs)
% 
% % specifying default values
% I.key1 = {1 2 3};
% I.key2 = 1:3;
% I.key3 = 'abc';
% [I, C] = parse_optInputs_keyvalue({'key1', {4,5,6}}, I)
% 
% % empty lists mean unspecified
% I.key1 = {1 2 3};
% I.key2 = 1:3;
% I.key3 = 'abc';
% [I, C] = parse_optInputs_keyvalue({'key1', {}}, I, 'empty_means_unspecified', true)
% 
% % use defaults to catch a spelling mistake
% I = parse_optInputs_keyvalue({'keys1', {4,5,6}}, I)

% 2016-08-27: Created, Sam NH
% 
% 2018-05-23: Added functionality to detect if an optional input has
% changed. Added functionality to have empty lists mean unspecified values.
% 
% 2019-01-18: Added C_value, useful for creating strings of parameters that
% have changed

clear P;
P.empty_means_unspecified = false;
P.ignore_bad_keys = false;
P.ignore_mismatch_class = {};
P.always_include = {};
P.always_exclude = {};
P.maxlen = 100;
P.delimiter = '/';
P.noloop = false;
P.paramstring = false;
if ~isempty(varargin)
    if strcmp(varargin{1}, 'noloop')
        P.noloop = true;
    else
        P = parse_optInputs_keyvalue(varargin, P, 'noloop');
        if isempty(P.empty_means_unspecified)
            P.empty_means_unspecified = false;
        end
    end
end

% should be an event number of arguments
n_optargs = length(optargs);
if mod(n_optargs,2)~=0
    error('There are not an even number of optional inputs');
end

% initialize with empty structure if not specified
if nargin < 2
    I = struct;
else % extract list of possible keys from the values of I if specified
    possible_keys = fieldnames(I);
end

% check keys are not repeated
try
    if length(unique(optargs(1:2:n_optargs))) ~= length(optargs(1:2:n_optargs))
        error('Duplicate keys');
    end
catch
    error('where am i');
    keyboard
end

% assume keys are not changed unless they are
if exist('possible_keys', 'var')
    for i = 1:length(possible_keys)
        C.(possible_keys{i}) = false;
    end
end

C_value = struct;

% immediately return if there are no optional arguments
if n_optargs == 0 && ~P.paramstring
    all_keys = {};
    return;
end

% assign keys and values
i_key = 1:2:n_optargs;
i_val = 2:2:n_optargs;
n_pairs = checkint(n_optargs/2);
all_keys = cell(1, n_pairs);
for j = 1:n_pairs
    
    key = optargs{i_key(j)};
    value = optargs{i_val(j)};
    all_keys{j} = key;

    % check key is a string
    if ~ischar(key)
        error('Optional arguments not formatted propertly\nAll keys must be strings\n');
    end
    
    % check key is one of the possible keys, and if it has the same type class
    % type
    if exist('possible_keys', 'var')
        if ~any(strcmp(key, possible_keys))
            if P.ignore_bad_keys
                continue;
            else
                error(['Optional arguments not formatted propertly\n' ...
                    '''%s'' not a valid key\n'], key);
            end
        end
        
        if ~isequal(class(I.(key)), class(value)) && ~any(ismember(P.ignore_mismatch_class, key))
            allowed_class_swaps = {'double', 'int32'};
            allowed_swap = false;
            for k = 1:size(allowed_class_swaps,1)
                if P.empty_means_unspecified && isempty(value)
                    allowed_swap = true;
                end
                if strcmp(class(I.(key)), allowed_class_swaps{k,1}) ...
                        && strcmp(class(value), allowed_class_swaps{k,2})
                    allowed_swap = true;
                end
                if strcmp(class(value), allowed_class_swaps{k,1}) ...
                        && strcmp(class(I.(key)), allowed_class_swaps{k,2})
                    allowed_swap = true;
                end                
            end
            if ~allowed_swap
                error(['Optional arguments not formatted propertly\n' ...
                    'Value of ''%s'' should be of type %s\n'], key, class(I.(key)));
            end
        end
    end
        
    % assign
    if P.empty_means_unspecified && isempty(value)
        % do nothing
    else
        if ~isfield(I, key) || ~isequaln(I.(key), value)
            I.(key) = value;
            C.(key) = true;
            C_value.(key) = I.(key);
        end
    end
end

% reorder fiels to match original structure
C_value = orderfields(C_value, intersect(fieldnames(I), fieldnames(C_value), 'stable'));

if P.paramstring && ~P.noloop
    paramstring = optInputs_to_string(I, C_value, P.always_include, P.always_exclude, ...
        'maxlen', P.maxlen, 'delimiter', P.delimiter);
else
    paramstring = '';
end
