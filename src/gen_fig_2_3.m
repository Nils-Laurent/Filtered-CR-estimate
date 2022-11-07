clear workspace;
close all;
 
cas = 1;
prefix = 'LC';
% gen_fig_q_filt;
test_LCR;
% return;

cas = 2;
prefix = 'poly';
% gen_fig_q_filt;
test_LCR;
% return;

cas = 3;
prefix = 'sin';
% gen_fig_q_filt;
test_LCR;

return;
cas = 1;
prefix = 'LC';
test_LCR;

cas = 2;
prefix = 'poly';
test_LCR;

cas = 3;
prefix = 'sin';
test_LCR;