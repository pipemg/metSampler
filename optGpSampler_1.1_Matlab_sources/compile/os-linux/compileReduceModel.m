function compileReduceModel()
    % This file compiles the C++ source code to a MATLAB .mexa file on the
    % LINUX operating system.In general this is only needed when the precompiled libraries do not
    % work.
    %
    % Requirements: 
    %  - A C++ compiler that support openMP (e.g. GCC or MS Visual Studio Pro/Ultimate)
    %  - Armadillo 4.200 or later
    %  - Boost C++ 1.55.0 or later 
    %  - One (or more) of the following solvers:
    %  1) IBM ILOG Cplex version 12.6 or later
    %  2) Gurobi version 5.6 or later
    %  3) GNU GLPK version 4.53 or later
    % 
    % NB: mind the starting space for the file and directory names!
    clc;

    
    % Source files
    cppFiles = [    ' ../src/CbModel.cpp', ...
                    ' ../src/CbModelCreator.cpp', ...
                    ' ../src/CbModelIOMatlab.cpp', ...
                    ' ../src/mxValidate.cpp', ...
                    ' ../src/reduceModel.cpp', ...
                    ' ../src/reduceModelMex.cpp'];                    

    cppIncludeDirs =   ' -I"../include"';
    outputFile = ' cReduceModel';


    %% Solver specific includes (at least 1 is needed)
    cplexIncludeDirs     =   [   ' -I"/home/wout/solvers/CPLEX_Enterprise_Server126/CPLEX_Studio/cpoptimizer/include"', ...
                                 ' -I"/home/wout/solvers/CPLEX_Enterprise_Server126/CPLEX_Studio/cplex/include"', ...
                                 ' -I"/home/wout/solvers/CPLEX_Enterprise_Server126/CPLEX_Studio/concert/include"'];

    % Specify path to cplex include directory, or if you do not want to link to IBM ILOG Cplex, uncomment the line below
    % cplexIncludeDirs = '';
    
    % Specify path to gurobi include directory, or if you do not want to link to GUROBI, uncomment the line below
    gurobiIncludeDirs    =  [    ' -I"/home/wout/solvers/gurobi562/linux64/include"'];
    % gurobiIncludeDirs = '';

    % Specify path to GLPK include dir, or if you do not want to link to GLPK, uncomment the line below
    glpkIncludeDirs      =      ' -I"/home/wout/solvers/glpk-4.53/lib/include"';
    % glpkIncludeDirs = '';

    
    %% Other mandatory includes
    % Path to ARMADILLO include
    armadilloIncludeDir         =       ' -I"/home/wout/software/armadillo-4.200.0/include"';
    
    % Path to Boost C++ include
    boostIncludeDir             =       ' -I"/home/wout/software/boost_1_55_0"';


    matlabLibs = fullfile(matlabroot, 'bin',  computer('arch'));
  
    cmd = ['mex -v -largeArrayDims COMPFLAGS="\$COMPFLAGS" CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" ', ...
            cppIncludeDirs, cppFiles, ' -output ', outputFile,  ...
            cplexIncludeDirs, gurobiIncludeDirs, glpkIncludeDirs, ... 
            armadilloIncludeDir, boostIncludeDir, ...
            ' -L"', matlabLibs, '" -lut -ldl'];
        
    eval(cmd);
    fprintf('reduceModel.%s successfully compiled! \n', mexext);
    system(sprintf('mv cReduceModel.%s ../../', mexext));
   
    fprintf('Files moved to package:: optGpSampler \n');
end




  

