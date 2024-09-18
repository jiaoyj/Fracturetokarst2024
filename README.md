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

7. Input files instruction:

   
             //////////////////instruction of input file: In_DFractureParabenchmark.dat as an example with input data
            2 1 3 // nboun1 nboun2 nbounf, which are the number of the first boundaries, the number of the second boundaries,the number of the third boundaries
            
            0 0 0 200        //xy12[i][1] xy12[i][2] xy12[i][3] xy12[i][4], one of the first boudaries, the coordinate of its two points(x,y unit m)
            
            2000 0 2000 200  //xy12[i][1] xy12[i][2] xy12[i][3] xy12[i][4], one of the first boudaries, the coordinate of its two points
            
            0  0 2000 0       //xy12[i][1] xy12[i][2] xy12[i][3] xy12[i][4], one of the second boudaries, the coordinate of its two points
            
            0 200 0 500         //xy12[i][1] xy12[i][2] xy12[i][3] xy12[i][4], one of the third boudaries, the coordinate of its two points
            
            2000 200 2000 500   //xy12[i][1] xy12[i][2] xy12[i][3] xy12[i][4], one of the third boudaries, the coordinate of its two points
            
            0 500 2000 500      //xy12[i][1] xy12[i][2] xy12[i][3] xy12[i][4], one of the third boudaries, the coordinate of its two points
            
            0 0           // the coordinate of each point on the boundaries, clockwise
            
            0 200         //
            
            0 500         //
            
            2000 500      //
            
            2000 200      //
            
            2000 0        //the coordinate of each point on the boundaries, clockwise
            
            200 200       // the water head of 2 first given water head boundaries
            
             0            // the flux of only one sencond  no flow boundary
            
            0 // the number of the determined fractures. When it is 0, there are only random fractures. The fracture network is composed of determined and random fractures. 
            
            -1 -44 201 6 0.005     //the coordinate of two points of a fracture and its aperture(cm)
            
            -1 -39 201 11 0.005     //
            
            -1 -34 201 16 0.005     //
            
            -1 -29 201 21 0.005     //
            
            ......                    //the remaining fractures
            
            
            
            /////////////////////instruction of input file: In_RFracturePara.dat ,which generate the random fractures
            
            0 2000 0 500  // minx,maxx,miny,maxy,the range of generating fracures
            
            20 5          //ndx,ndy, the number of subareas in x and y directions
            
            2             // the number of groups of generated random fracutres
            
            3 240 5 220 260               //ndisl[i], meanl[i], sigmal[i], minl[i], maxl[i], the fracture length parameters of the first group of fractures, distribution, mean, standard deviation, minimum, maximum
            
            4 0.01 0.001 0.008 0.012      //ndisb[i], meanb[i], sigmab[i], minb[i], maxb[i], the aperture parameters:distribution, mean, standard deviation, minimum, maximum
            
            2 11 1 9 13                   //ndisd[i], meand[i], sigmad[i], mind[i], maxd[i], the direction parameters:distribution, mean, standard deviation, minimum, maximum
            
            20                            //the parameter of fracute density of the first group
            
            3 240 5 220 260               //ndisl[i], meanl[i], sigmal[i], minl[i], maxl[i],the fracture parameters of the second group of fractures
            
            4 0.01 0.001 0.008 0.012      //ndisb[i], meanb[i], sigmab[i], minb[i], maxb[i]
            
            2 85 2 82 88                  //ndisd[i], meand[i], sigmad[i], mind[i], maxd[i]
            
            10                            //the parameter of fracute density of the second group
            
            
            ///////////////////////instruction of input file: In_TimeStep.dat
            
            100  1 500 // the number of time steps, output control, and a time step(year)
            
            
            ////////////////////////instruction of input file: In_InfilPara.dat
            
            1000.0  0.2  // The annual rainfall recharge (mm) and the coefficient of the rainfall infiltration
            
            
            







