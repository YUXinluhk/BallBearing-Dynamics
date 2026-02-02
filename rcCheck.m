function res = rcCheck(Q_o, a_o, E_o, trig, Q_i, a_i, E_i)
Q_o_sort=sort(Q_o);
rep_value=Q_o_sort(round(length(Q_o)/2));
rep_element=find(Q_o==rep_value);
rep_element=rep_element(1);
Q_o_med=Q_o(rep_element);
Q_i_med=Q_i(rep_element);
a_o_med=a_o(rep_element);
a_i_med=a_i(rep_element);
E_i_med=E_i(rep_element);
E_o_med=E_o(rep_element);
trig_med=trig(rep_element);
res=(Q_o_med*a_o_med*E_o_med*trig_med<=Q_i_med*a_i_med*E_i_med);
end
