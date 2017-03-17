function vars = gpu_whos()
% vars = gpu_whos()
%
% An enhanced whos to list gpuArray variables and AOSim2 objects that use
% GPU resources.

vars_ = evalin('base','whos()');
fprintf('Found %d base variables.\n',length(vars_));

nn = 1;
for n=1:length(vars_)
    
    uses_gpu = false;
    if(strcmp(vars_(n).class,'gpuArray'))
        uses_gpu = true;
    end
    
    if(strncmp(vars_(n).class,'AO',2))
        uses_gpu = evalin('base',sprintf('%s.useGPU',vars_(n).name));
    end
    
    if(uses_gpu)
        fprintf('%s\t\t%s ',...
            vars_(n).name, vars_(n).class);
        
%         if(strcmp(vars_(n).class,'gpuArray'))
%             fprintf('(GPU) ')
%         end
        
%         if(strncmp(vars_(n).class,'AO',2))
%             fprintf('(AOSim2) ')
%         end
        fprintf('\n')
        
        %if(strcmp(vars_(n).class,'gpuArray') || strncmp(vars_(n).class,'AO',2))
        vars(nn) = vars_(n);
        nn = nn + 1;
    end
end

