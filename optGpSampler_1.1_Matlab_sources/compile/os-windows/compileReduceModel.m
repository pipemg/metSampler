function compileReduceModel()
    % This file compiles the C++ source code to a MATLAB .mexa file on the
    % WINDOWS operating system.In general this is only needed when the precompiled libraries do not
    % work.
    %
    % Requirements: 
    %  - A C++ compiler that support openMP (e.g. GCC or MS Visual Studio Pro\Ultimate)
    %  - Armadillo 4.200 or later
    %  - Boost C++ 1.55.0 or later 
    %  - One (or more) of the following solvers:
    %  1) IBM ILOG Cplex version 12.6 or later
    %  2) Gurobi version 5.6 or later
    %  3) GNU GLPK version 4.53 or later
    % 
    % NB: mind the starting space for the file and directory names!
    
    clear mex;
    clc;

    % Source files
    cppFiles = [    ' ..\src\CbModel.cpp', ...
                    ' ..\src\CbModelCreator.cpp', ...
                    ' ..\src\CbModelIOMatlab.cpp', ...
                    ' ..\src\mxValidate.cpp', ...
                    ' ..\src\reduceModel.cpp', ...
                    ' ..\src\reduceModelMex.cpp'];                    

    cppIncludeDirs =   ' -I"..\include"';
    outputFile = ' cReduceModel';


    %% Solver specific includes (at least 1 is needed)
    cplexIncludeDirs     =   [  ' -I"D:\solvers\CPLEX_Enterprise_Server126\CPLEX_Studio\cpoptimizer\include"', ...
                                ' -I"D:\solvers\CPLEX_Enterprise_Server126\CPLEX_Studio\cplex\include"', ...
                                ' -I"D:\solvers\CPLEX_Enterprise_Server126\CPLEX_Studio\concert\include"'];

    % Specify path to cplex include directory, or if you do not want to link to IBM ILOG Cplex, uncomment the line below
    % cplexIncludeDirs = '';
    
    % Specify path to gurobi include directory, or if you do not want to link to GUROBI, uncomment the line below
    gurobiIncludeDirs    =  [    ' -I"D:\solvers\gurobi560\win64\include"'];
    % gurobiIncludeDirs = '';

    % Specify path to GLPK include dir, or if you do not want to link to GLPK, uncomment the line below
    glpkIncludeDirs      =      ' -I"D:\Dropbox\software\glpk-4.53\src"';
    % glpkIncludeDirs = '';

    
    %% Other mandatory includes
    % Path to ARMADILLO include
    armadilloIncludeDir         =       ' -I"D:\Dropbox\software\armadillo-4.200.0\include"';
    
    % Path to Boost C++ include
    boostIncludeDir             =       ' -I"D:\Dropbox\software\boost_1_55_0"';
    matlabLibs = fullfile(matlabroot, 'bin',  computer('arch'));
  
    cmd = ['mex -v -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp"  CXXFLAGS="$CXXFLAGS /openmp" LDFLAGS="$LDFLAGS /openmp" ', ...
            ' LINKFLAGS="$LINKFLAGS /DELAYLOAD:CbModelGurobi.dll /DELAY:UNLOAD /DELAYLOAD:CbModelCplex.dll /DELAY:UNLOAD /DELAYLOAD:CbModelGlpk.dll /DELAY:UNLOAD  "',  ...            
            boostIncludeDir, armadilloIncludeDir, cplexIncludeDirs, gurobiIncludeDirs, glpkIncludeDirs, ... 
            cppIncludeDirs, cppFiles, ' -output ', outputFile, ...
            ' -L"..\..\windows_lib" -L"', matlabLibs, '" -lCbModelCplex -lCbModelGurobi -lCbModelGlpk -lmwblas -lmwlapack -lut -ldelayimp'];
        
    eval(cmd);
    fprintf('reduceModel.%s successfully compiled! \n', mexext);
    system(sprintf('move cReduceModel.%s ..\\..', mexext));
   
    fprintf('Files moved to package:: optGpSampler \n');
end