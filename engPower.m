function res = engPower(b, e)
%This is a hack that vastly improves the convergence of the model:
%Since the force model is F~delta^(3/2), if the solver ever guesses a
%negative delta, it will crash. With this hack, it will be led back to
%a correct solution. CAUTION, however: This class cannot desctibe cases
%where one or more balls are truely out of contact. Those cases will be
%indicated by a negative delta, but the results produced here are
%useless.
res = sign(b).*(abs(b).^e);
end
