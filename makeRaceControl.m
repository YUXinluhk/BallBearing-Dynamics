function [rc_id, lambda_i, lambda_o]=makeRaceControl(scenario)
switch scenario
    case 'inner'
        rc_id=0;
        lambda_i=1;
        lambda_o=1;
    case 'outer'
        rc_id=1;
        lambda_i=0;
        lambda_o=2;
    otherwise
        error(['Invalid Race control mode.'...
            ' Valid modes: inner, outer'])
end
end
