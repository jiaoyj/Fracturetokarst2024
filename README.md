# Fracturetokarst2024

   Karst fracture flow and evolution under rainfall conditions

1. Fracturetokarst program is originally developed by Professor YU Qingchun, which is a C++ code. the original program on Windows has friendly GUI interface, which can be free available on the Rock Fracture Group website http://www.rockfractures.com/. 
2. Recently, considering the special circumstances of fracture phreatic water table and rainfall recharge, it was improved and compiled under Linux 64bit environment and run at high performance computing platform. Since the new developed version is under Linux and the running time has increased significantly, the exe on Windows is not suitable and the friendly GUI interfacehas not developed with the latest code.

   It need a linear algebra libary Armadillo as C++ libary to do matrix computing. We should install Armadillo package on scientific computing platform first before compiling the source code. The download website of Armadillo package is https://arma.sourceforge.net/.

3. The source code files are main.cpp FractureModel.cpp  KarstModel.cpp  RandomGenerator.cpp and their .h files. So for example, when solving the free surface and karst evolution in the paper, the compiling source code can be as following: g++ main.cpp RandomGenerator.cpp FractureModel.cpp KarstModel.cpp -o karstprofile -L /public/home/jiaoyj/software/armadillo9/lib64 -larmadillo.

   So the exe karstprofile  is mainly aimed to solving the free surface and its evolution under rainfall conditions.

4. The file run.slurm is the job file asigning computing nodes and cores, which start the karstprofile exe and inputfiles. the command line on the scientific computing platform could be as：sbatch run.slurm. It depends on the computing platforms.

5. The Input files are In_DFractureParabenchmark.dat  In_InfilPara.dat  In_RFracturePara.dat  In_TimeStep.dat.

6. The output files will be created with out_ prefix in the file name.



