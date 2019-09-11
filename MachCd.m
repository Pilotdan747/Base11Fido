clear
clc
close all

cd = 0.0087;
cdmm = 0.51;
cdm = 0.3;

maxMach = 4;
step = 0.01;
iMax = maxMach/step;

Cd = zeros(1, iMax);

for i = 1:iMax
    if (i < 0.8/step)
        Cd(i) = cd;
    end
    
    if (i >= 0.8/step && i < 1/step)
        Cd(i) = 1.0273*(i*step)^3 - 0.5173; %numbers from polyfit([0.8 1],[0.0087 0.51], 3)
    end
    
    if (i == 1/step)
        Cd(i) = cdmm;
    end
    
    if (i > 1/step && i <= 1.2/step)
        Cd(i) =  -0.2885*(i*step)^3 + 0.7985;
    end
    
    if (i > 1.2/step)
       Cd(i) = cdm;   
    end  
end

M = 0:step:maxMach-step;

plot(M, Cd)