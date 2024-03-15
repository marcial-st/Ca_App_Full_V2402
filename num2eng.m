function varargout = num2eng(num,varargin)
%NUM2ENG converts numbers to engineering notation.
%   Engineering notation is similar to scientific notation, except that
%   the exponent power of 10 is always an integer multiple of 3, and the
%   mantissa is scaled so that it lies inside [0,1000). The range of
%   exponents for which a SI-prefix is available goes from 10^-24 to 10^24.
%   
%   Examples (of scalar inputs):
%       - num2eng(4700)
%         ans = '4.7 k'
%       
%       - num2eng(82000,[1 4], 'FullName',true)
%         ans = '82 kilo'
%       
%       - num2eng(1000,4)
%         ans = '1.000 k'
%       
%       - num2eng(-3.527e6,3)
%         ans = '-3.53 M'
%       
%       - num2eng(9999,[1 2])
%         ans = '10 k'
%       
%       - num2eng(999,3)
%         ans = '999'
%       
%       - num2eng(100, 'ForceExponent',3)
%         ans = '0.1 k'
%       
%       - [ 'Distance: ', num2eng(4.2e4, 'FullName',true, 'Unit','meter') ]
%         ans = 'Distance: 42 kilometers'
%       
%       - [ 'Distance: ', num2eng(1e3, 'FullName',true, 'Unit','meter') ]
%         ans = 'Distance: 1 kilometer'
%   
%   Syntax:
%       eng = num2eng(num) returns a string containing the engineering
%       notation of 'num'. 'num' can be any numerical scalar, vector or
%       matrix. If the input isn't scalar the output is a cellstring with
%       the same shape as in the input. The number of significant figures
%       defaults to up to 5 ( 'sigfigs'=[1 5] ). The output contains a
%       space between the number and the prefix part, but doesn't have
%       trailing whitespaces.
%       
%       __ = num2eng(num,sigfigs) gives control over the number of
%       significant figures in the output. sigfigs must be positive. It can
%       either be a scalar or two-element vector describing [min max]. In
%       the latter case trailing zeros are, if applicable, omitted.
%       
%       __ = num2eng(__,Name,Value,...) allows the user to specify one or
%       more of these options to take effect:
%         
%         - 'SynchronizeExponents' : (false) / 'min' / 'max' / 'median'
%           Causes the outputs to all have the same SI-prefix/exponent. The
%           value for this ('min','max' or 'median') determines what
%           SI-prefix is chosen when values in num are of varying orders of
%           magnitude.
%         
%         - 'ForceExponent' : (false) / [exp] / [min max]
%           Forces all outputs to have a given exponent or confines them
%           within a given minimum and maximum exponent range. Exponent
%           values must be multiples of 3. If this option is present,
%           'SynchronizeExponents' is ignored.
%         
%         - 'DistinctOutput' : (false) / true
%           Ensures distinct outputs. The maximum number of significant
%           figures is increased as needed to achieve this.
%         
%         - 'OmitSpace' : (false) / true
%           Omits the whitespace between the number and the SI-prefix. This
%           option is ignored if there are two output arguments.
%         
%         - 'FullName' : (false) / true
%           Writes prefixes as full names instead of one-letter
%           abbreviations.
%         
%         - 'Unit' : (false) / 'unit' / {'unit1','unit2',...}
%           Appends the given unit string(s) to the output. If a cell-array
%           is provided, it must have the same number of elements as the
%           'num' input. If necessary 'reshape' is called to match the
%           shape of 'num'. Each unit string gets processed by 'strtrim'.
%           Written out units should be given in singular form. Whenever
%           the printed value is not '1', a plural 's' is appended to the
%           unit string. If there are two output arguments and the exponent
%           is forced to a single value ( via the option 'ForceExponent'
%           with a single exponent or via 'SynchronizeExponents' ), '(s)'
%           or 's' is appended to the second output, depending on if some,
%           but not all, or all printed values are not '1'. No plural 's'
%           is appended if the unit string already ends with a 's'.
%         
%         - 'NoPlural' : (false) / true
%           Disables the automatic appending of a plural 's' to the unit.
%           This option is ignored if full names aren't effective aswell.
%         
%         All options default to false, which disables them. The option
%         names can be abbreviated as long as they remain unambigious.
%         A working set of short names is 'sync', 'force', 'distinct',
%         'omit', 'full', 'unit' and 'nop'. The shortest possible set of
%         names is 's', 'fo', d', 'fu', 'u' and 'n'.
%       
%       [num,sip] = num2eng(num,__) returns number and SI-prefix part
%       separately. Neither part contains a whitespace. If the input 'num'
%       isn't scalar, the output 'num' is a cellstring. 'sip' is a
%       cellstring if the input 'num' isn't scalar and the exponent in the
%       output isn't forced to a single value ( via the options
%       'SynchronizeExponents' or 'ForceExponent',[exp] ).
%   
%   Performance Considerations:
%       - The complexity of this function means it has considerable
%         overhead for simple tasks and especially scalar inputs. To help
%         mitigate this use this function in a vectorized manner where
%         possible. The effect of the overhead diminishes for big input
%         arrays.
%       
%       - Notice that the 'DistinctOutput' option can be especially costly
%         if multiple iterations of increasing the number of significant
%         figures are required to achieve distinctness. Choose a generous
%         upper bound for the number of significant figures if you expect
%         input values to be close to each other.
%   
%   Miscellaneous notes / side effects:
%       - 'warning' modes get reset to their default state
%         ( 'backtrace'->'on', 'verbose'->'off' )
%   
%   Version 1.2.4.3
%   Created by Roman Mueller-Hainbach, 29th of November 2017
%   Based on work by Rick Rosson and Stephen Cobeldick

%% Validate inputs

% Validate number of outputs
nargoutchk(0,2)

% Declare peristent parser variable
persistent parser

% Validate 'num' argument
validateattributes( num, {'numeric'}, {'real','nonsparse'}, ...
                    mfilename, 'Number', 1 )

if ~isempty(varargin) && isnumeric(varargin{1})
    % Retrieve 'SignificantFigures' argument
    sigfigs = varargin{1};
    % Validate and prepare 'sigfigs' parameter
    validateattributes( sigfigs, {'numeric'}, ...
                        { 'nonempty', 'positive', 'vector', ...
                          'finite', 'real', 'increasing', ...
                          'nonsparse' }, mfilename, ...
                        'SignificantFigures', 2 )
    assert( numel(sigfigs) <= 2, [ '''SignificantFigures'' must be a ' ...
                                   'scalar or two-element vector.' ] )
    if isscalar(sigfigs)
        sigfigs(2) = sigfigs;
    end
    varargin(1) = [];
else
    sigfigs = [1 5];
end

% Parse options
if ~isempty(varargin)
    
    % Create parser object
    if isempty(parser)
        
        parser = inputParser;
        parser.FunctionName = mfilename;
        
        % Define default value for parameters
        PARAM_DEF = false;
        
        % Add 'SynchronizeExponents' parameter
        syncValidFcn = @(s)validateattributes( s, {'logical','char'}, ...
                                               {'vector','nonempty'}, ...
                                               mfilename, ...
                                               'SynchronizeExponents' );
        parser.addParameter( 'SynchronizeExponents', PARAM_DEF, ...
                             syncValidFcn )
        
        % Add 'DistinctOutput' parameter
        distValidFcn = @(dis)validateattributes( dis, {'logical'}, ...
                                                 {'scalar','nonempty'}, ...
                                                 mfilename, ...
                                                 'DistinctOutput' );
        parser.addParameter( 'DistinctOutput', PARAM_DEF, distValidFcn )
        
        % Add 'OmitSpace' parameter
        omitValidFcn = @(omi)validateattributes( omi, {'logical'}, ...
                                                 {'scalar','nonempty'}, ...
                                                 mfilename, 'OmitSpace' );
        parser.addParameter( 'OmitSpace', PARAM_DEF, omitValidFcn )
        
        % Add 'FullName' parameter
        fullValidFcn = @(ful)validateattributes( ful, {'logical'}, ...
                                                 {'scalar','nonempty'}, ...
                                                 mfilename, 'FullName' );
        parser.addParameter( 'FullName', PARAM_DEF, fullValidFcn )
        
        % Add 'NoPlural' parameter
        nopValidFcn = @(ful)validateattributes( nop, {'logical'}, ...
                                                 {'scalar','nonempty'}, ...
                                                 mfilename, 'NoPlural' );
        parser.addParameter( 'NoPlural', PARAM_DEF, nopValidFcn )
        
        % Add 'ForceExponent' parameter
        forceValidFcn = @(forc)validateattributes( forc, { 'logical', ...
                                                           'numeric' }, ...
                                                   { 'nonempty', ...
                                                     'vector'  }, ...
                                                   mfilename, ...
                                                   'ForceExponent' );
        parser.addParameter( 'ForceExponent', PARAM_DEF, forceValidFcn )
        
        % Add 'Unit' parameter
        unitValidFcn = @(u)validateattributes( u, { 'logical', 'char', ...
                                                    'cell' }, ...
                                                  { 'nonempty', ...
                                                    'nonsparse' }, ...
                                                  mfilename, 'Unit' );
        parser.addParameter( 'Unit', PARAM_DEF, unitValidFcn )

    end % if
    
    parser.parse(varargin{:})
    parsed = parser.Results;
    
    % Retrieve 'SynchronizeExponents' argument from parser
    sync = parsed.SynchronizeExponents;
    % Validate and prepare 'sync' parameter
    if islogical(sync)
        assert( isequal(sync,false), ...
                [ 'Expected ''SynchronizeExponents''' ...
                  'to be false or a string.' ] )
    else % ischar(sync)
        sync = validatestring( sync, {'min','max','median'}, ...
                               mfilename, 'SynchronizeExponents' );
    end
    
    % Retrieve 'DistinctOutput' argument from parser
    distinct = parsed.DistinctOutput;
    % Validate and prepare 'shouldBeDistinct' parameter
    if distinct
        isDistinct = ~any( diff( sort(num) ) == 0 );
        assert( isDistinct, 'Distinct output requires a distinct input.' )
    end
    
    % Retrieve 'FullName' argument from parser
    full = parsed.FullName;
    
    % Modify warning modes
    warning('off','backtrace')
    warning('off','verbose')
    
    try
    
        % Retrieve 'NoPlural' argument from parser
        noplural = parsed.NoPlural;
        if noplural && ~full
            warning( 'NUM2ENG:ignored_option:noplural', ...
                     'Option ''NoPlural'' ignored.' )
        end
        
        % Retrieve 'ForceExponent' argument from parser
        force = parsed.ForceExponent;
        % Validate and prepare 'force' parameter
        if islogical(force)
            assert( isequal(force,false), ...
                    [ 'Expected ''ForceExponent'' ' ...
                      'to be false or numerical.' ] )
        else % isnumeric(force)
            validateattributes( force, {'numeric'},                  ...
                                { 'real','increasing','nonsparse' }, ...
                                mfilename, 'ForceExponent'           )
            assert( numel(force) <= 2, ...
                    [ '''ForceExponent'' must be a scalar ' ...
                      'or two-element vector.' ] )
            assert( all( mod(force,3) == 0), ...
                    'Forced exponent(s) must be a multiple of 3.' )
            % Check applicability of parameter 'sync'
            if ~isequal(sync,false)
                warning( 'NUM2ENG:ignored_option:sync', ...
                         'Option ''SynchronizeExponents'' ignored.' )
                sync = false;
            end
        end

        % Retrieve 'OmitSpace' argument from parser
        omit = parsed.OmitSpace;
        % Check applicability of parameter 'omit'
        if nargout == 2 && omit
            warning( 'NUM2ENG:ignored_option:omit', ...
                     'Option ''OmitSpace'' ignored.' )
        end
        
    catch exception
        % Restore default warning modes before exiting function
        % by throwing an error
        warning('on','backtrace')
        warning('off','verbose')
        rethrow(exception)
    end
    
    % Restore default warning modes
    warning('on','backtrace')
    warning('off','verbose')
    
    % Retrieve 'Unit' argument from parser
    unit = parsed.Unit;
    % Validate and prepare 'unit' parameter
    if islogical(unit)
        assert( unit == false, ...
                'Expected ''Unit'' to be false or a (cell)string.' )
    else
        if iscell(unit)
            assert( numel(unit) == numel(num), ...
                    [ 'The number of ''unit'' strings doesn''t ' ...
                      'match the number of inputs.' ] )
            unit = reshape( unit, size(num) );
        end
        unit = strtrim(unit);
    end
    
else
    
    % Assign default values to unparsed parameters
    sync = false;
    distinct = false;
    full = false;
    force = false;
    omit = false;
    unit = false;
    
end % if

%% Preallocate output
if isscalar(num)
    varargout{1} = {};
    if nargout == 2
        varargout{2} = {};
    end
    exponentIsVariable = false;
else
    varargout{1} = cell(size(num));
    exponentIsVariable = ~isscalar(force) || ...
                         ( islogical(force) && isequal(sync,false) );
    if nargout == 2
        if exponentIsVariable
            varargout{2} = cell(size(num));
        else
            varargout{2} = '';
        end
    end
end

%% Preprocess num
% Determine signs of num and remove them from the numbers themselves
isnegative = num < 0;
num = abs(num);
% Round numbers with the maximum number of significant digits
num = round( num, sigfigs(end), 'significant' );

% Determine exponents
exponents = 3 .* floor( floor(log10(num)) ./ 3 );
exponents( exponents == -Inf ) = 0;

% Manipulate exponents according to passed options
if isnumeric(force)
    if isscalar(force)
        [exponents(:)] = force;
    else
        exponents = min( max(exponents,force(1)), force(2) );
    end
elseif ~isequal(sync,false)
    value = feval(sync,exponents);
    [exponents(:)] = value - rem(value,3);
end

%% Convert each number
for i = 1:numel(num)
    
    number = num(i);
    
    % Look at sign of number
    if isnegative(i)
        signstr = '-';
    else
        signstr = '';
    end
    
    % Handle Inf and NaN values
    if isinf(number) || isnan(number)
        if isnan(number)
            str = 'NaN';
        else % isinf
            str = 'Inf';
        end
        varargout{1}{i} = [ signstr, str ];
        if nargout == 2 && exponentIsVariable
            varargout{2}{i} = '';
        end
        continue
    end
    
    % Determine mantissa
    mantissa = number / 10^exponents(i);
    
    % Determine number of digits in whole number and fractional part
    if mantissa == 0.0
        mantissaoom = 0;
    else
        mantissaoom = floor( log10(mantissa) );
    end
    wholedigs = max( 0, mantissaoom + 1 );
    if sigfigs(2) ~= sigfigs(1)
        int = uint64( mantissa * 10^sigfigs(2) );
        if int == 0
            fracdigs = sigfigs(1)-1;
        else
            reqdigs = max([0,find( sprintf('%i',int) ~= '0', 1, 'last' )]);
            if wholedigs > 0
                fracdigs = max( 0, max(sigfigs(1),reqdigs) - wholedigs );
            else
                fracdigs = -mantissaoom;
            end
        end
    else
        fracdigs = max(0, sigfigs(2) - wholedigs );
    end
    
    % Determine SI-prefix for exponent
    expOutOfRange = false;
    if mantissa == 0 && isequal(sync,false) && islogical(force)
        expstr = '';
    elseif full
        switch exponents(i)
            % Ordered by commonness
            case   0; expstr =  '';
            case   3; expstr = 'kilo';
            case  -3; expstr = 'milli';
            case   6; expstr = 'mega';
            case  -6; expstr = 'micro';
            case  -9; expstr = 'nano';
            case -12; expstr = 'pico';
            case   9; expstr = 'giga';
            case -15; expstr = 'femto';
            case  12; expstr = 'tera';
            case  15; expstr = 'peta';
            case -18; expstr = 'atto';
            case  18; expstr = 'exa';
            case -21; expstr = 'zepto';
            case  21; expstr = 'zetta';
            case -24; expstr = 'yocto';
            case  24; expstr = 'yotta';
            otherwise
                expstr = sprintf( 'e%g ', exponents(i) );
                expOutOfRange = true;
        end % switch
    else
        switch exponents(i)
            % Ordered by commonness
            case   0; expstr =  '';
            case   3; expstr = 'k';
            case  -3; expstr = 'm';
            case   6; expstr = 'M';
            case  -6; expstr = 'u';
            case  -9; expstr = 'n';
            case -12; expstr = 'p';
            case   9; expstr = 'G';
            case -15; expstr = 'f';
            case  12; expstr = 'T';
            case  15; expstr = 'P';
            case -18; expstr = 'a';
            case  18; expstr = 'E';
            case -21; expstr = 'z';
            case  21; expstr = 'Z';
            case -24; expstr = 'y';
            case  24; expstr = 'Y';
            otherwise
                expstr = sprintf( 'e%g ', exponents(i) );
                expOutOfRange = true;
        end % switch
    end % if
    
    % Consider whitespace
    if expOutOfRange && ~isequal(unit,false)
        expstr(end) = [];
    end
    if omit || expOutOfRange || ( isempty(expstr) && isequal(unit,false) )
        spacestr = '';
    else
        spacestr = ' ';
    end
    
    % Assemble output string
    if nargout == 2
        varargout{1}{i} = [ signstr, ...
                            sprintf( '%0.*f', fracdigs, mantissa ) ];
        if exponentIsVariable
            varargout{2}{i} = expstr;
        end
    else
        varargout{1}{i} = [ signstr, ...
                            sprintf( '%0.*f', fracdigs, mantissa ), ...
                            spacestr, expstr ];
    end
    
    % Check for distinctness of strings
    if distinct
        eqs = strcmp( varargout{1}{i}, varargout{1}(1:i-1) );
        if any(eqs) && ( nargout <= 1 || ...
                         strcmp( varargout{2}(i), varargout{2}(eqs) ) )
            [varargout{:}] = feval( mfilename, parsed.Number, ...
                                    sigfigs+[0 1], 'distinct',true, ...
                                    'sync',sync, 'force',false, ...
                                    'unit',unit, 'nop',noplural, ...
                                    'omit',omit, 'full',full );
            return
        end
    end
    
end % for i

% Assign second output now if it is wasn't before
exponentIsFixed = ( ~islogical(force) && isscalar(force) ) || ...
                  ~isequal(sync,false);
if nargout == 2 && exponentIsFixed
    varargout{2} = expstr;
end

% Append unit string
if ~isequal(unit,false)
    iArgout = max(1,nargout);
    [varargout{iArgout}] = strcat( varargout{iArgout}, unit );
    % Append plural 's'
    iPlural = ~cellfun( @(str)str(end) == 's' || str(1) == '1' && ...
                        ( isscalar(str) || ...
                        ~any( str(2) == [46,48:57] ) ), varargout{1} );
    if ~noplural && any(iPlural)
        if iscell(varargout{iArgout})
            s = cell(size(varargout{iArgout}));
            s(iPlural) = {'s'};
            s(~iPlural) = {''};
            [varargout{iArgout}] = strcat( varargout{iArgout}, s );
        elseif all(iPlural)
            varargout{iArgout}(end+1) = 's';
        else
            varargout{iArgout} = [varargout{iArgout} '(s)'];
        end
    end
end

% Unwrap cell if output is scalar
if isscalar(varargout{1})
    varargout{1} = varargout{1}{1};
end
if nargout == 2 && iscell(varargout{2}) && isscalar(varargout{2})
    varargout{2} = varargout{2}{1};
end