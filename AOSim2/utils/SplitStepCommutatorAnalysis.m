function [Cn2,THICK,PHASE_RMS] = SplitStepCommutatorAnalysis(TOTAL_RANGE,R0_at_RANGE)
% [Cn2,THICK,PHASE_RMS] = SplitStepCommutatorAnalysis(TOTAL_RANGE,R0_at_RANGE)

PS = AOScreen(1);

Cn2 = PS.setR0(R0_at_RANGE,TOTAL_RANGE).Cn2;
THICK = 10.^(0:.1:ceil(log10(TOTAL_RANGE)));

FRIED = zeros(size(THICK));

for n=1:length(THICK)
    FRIED(n) = PS.setCn2(Cn2,THICK(n)).r0;
end

THETAs = PS.lambdaRef./FRIED;
PHASE_RMS = sqrt(6.88*(THETAs.*THICK./FRIED).^(5/3) .* (10e3 ./ THICK) );

clf;
% loglog(THICK,PHASE_RMS,'o-',THICK,1.1e-2*(THICK/100).^(4/3),'r--');
loglog(THICK,PHASE_RMS,'-','LineWidth',2);
grid;
biglabels('Split-Step (m)','rms phase uncertainty');
bigtitle(sprintf('Split-Step Commutator Change for r_0=%.1f cm @ %.0f km path',R0_at_RANGE*100,TOTAL_RANGE/1e3 ));
