function [valor,unidad] = eng2num(string_valor)
% [valor,unidad] = eng2num(string_valor)
% [0.01,'F'] = eng2num('10mF')
% [0.01,'']  = eng2num('10m')
% 0.01       = eng2num('10m')
%
% http://robotica.udl.es (05/2005)

% pasar a numero
% localizar el valor
for (i=1:1:length(string_valor))
    if (isempty(str2num(string_valor(1:i)))),
        i = i-1;
        break;
    end,
end,
valor = str2num(string_valor(1:i));

if (i < length(string_valor))
    i = i+1;
    % multiplicador
    switch string_valor(i:i),
    case 'P',
        valor = valor *1e15;
    case 'T',
        valor = valor *1e12;
    case 'G',
        valor = valor *1e9;
    case 'M',
        valor = valor *1e6;
    case 'K',
        valor = valor *1e3;
    case 'm',
        valor = valor *1e-3;
    case 'u',
        valor = valor *1e-6;
    case 'p',
        valor = valor *1e-9;
    case 'f',
        valor = valor *1e-12;
    otherwise,
        if (i > 1),
            i = i-1;
        end,
    end,
end,

% extraer unidad
unidad = [];
if (i < length(string_valor)),
    i = i+1;
    unidad = string_valor(i:end);
end,