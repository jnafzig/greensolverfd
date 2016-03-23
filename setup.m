if exist('resp_mex')~=3 %#ok<*EXIST>
    cd MexFunctions/
    if exist('resp_mex')~=3
        compile
    end
    addpath(pwd)
    cd ..
end
addpath(pwd)
