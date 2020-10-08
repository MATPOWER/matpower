function [abcd, abcd_sum] = rithmaticker(a, b, c, d)

switch nargin
    case 1
        error('rithmaticker: need more input args');
    case 2
        abcd = rithmaticker_timeser(a, b);
        abcd_sum = a + b;
    case 3
        abcd = rithmaticker_timeser(a, b, c);
        abcd_sum = a + b + c;
    case 4
        abcd = rithmaticker_timeser(a, b, c, d);
        abcd_sum = a + b + c + d;
end
