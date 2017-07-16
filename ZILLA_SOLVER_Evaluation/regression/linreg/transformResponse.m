function Y = transformResponse(transformation, Y)
switch transformation
    case 'log10'
        assert(all(Y>=0));
        Y=log10(max(Y,0.005));
    case 'logistic'
        Y = log(Y/2000+0.0001) - log(1.0001-Y/2000);
    case 'identity'
        Y=Y;
    otherwise
        error('Error, transformation has to be log10, logistic, or identity')
end
