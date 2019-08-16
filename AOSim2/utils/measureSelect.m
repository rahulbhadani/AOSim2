function [s,ang] = measureSelect()

% function [s,ang] = measureSelect()
% 
% 20081002: JLCodona

PNTS = ginput(2);

V = PNTS(2,:)-PNTS(1,:);

s = norm(V);
ang = atan2(V(2),V(1));

if(nargout==0) 
	fprintf('spacing=%.1f angle=%.1f deg\n', s,ang*180/pi);
end

return;
