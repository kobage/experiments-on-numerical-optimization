The CG_DESCENT-C-6.8 project is an adaptation of lcg for Visual Studio. The code is obtained from http://people.clas.ufl.edu/hager/files/CG_DESCENT-C-6.8.tar_.gz. We utilize this project to compare and evaluate our implemented solvers. Additionally, we have extended the project with a collection of test functions for unconstrained optimization, including a neural network model consisting of 10 neurons corresponding to the mnist dataset.

The MHB_01 project describes the experiments conducted in the paper titled "Speeding up the Convergence of the Polyak's Heavy Ball Algorithm" by Koba Gelashvili, Irina Khutsishvili, Luka Gorgadze, and Lela Alkhazishvili, published in Transactions of A. Razmadze Mathematical Institute, Volume 172 (2018), Issue 2, pages 176-188.

The L-BFGS_LB and MHB_LB projects describe the experiments conducted in the paper titled "Add-in for Solvers of Unconstrained Minimization to Eliminate Lower Bounds of Variables by Transformation" by Koba Gelashvili, published in Transactions of A. Razmadze Mathematical Institute, Volume 173 (2019), pages 39-46. The L-BFGS_LB project includes our C++ implementation of the L-BFGS algorithm. The line search algorithm is mostly consistent across all the projects.

The most interesting aspect of the OSC_FREE_MOMENTUMS project is the investigation of preconditioning in MHB. The experiments involve utilizing an incomplete Cholesky expansion with limited memory (C.J. Lin and J. J. Moré. Incomplete Cholesky factorizations with limited memory. SIAM J. on Scientific Computing, 21(1):24–45, 1999). The sparse Hessians are stored using std::vector<std::map<double>>, which simplifies coding but also introduces computational overhead. When applying the preconditioner to the 7 numerical test functions where lcg outperformed MHB in the experiments conducted in MHB_01 project, the results of preconditioned MHB significantly outperformed lcg for 6 of the functions. Only one case exhibited deterioration, accompanied by an increase in the number of iterations, indicating that the deterioration was not solely due to the sparse matrix implementation. These experiments demonstrate the potential of preconditioning in momentum-based optimizers, although achieving an efficient and robust implementation remains a challenge.

The Perf_profiles folder contains Python code that facilitates reading data from two or three files to construct performance profiles.

The RawCollectionOfMinimzers folder includes MHB (modified heavy ball algorithm), MNAG (Modified Nesterov Accelerated Gradient Algorithm), L-BFGS, ADAM, a line search algorithm, a collection of unconstrained minimization tests, and a logistic regression model. The project is developed in Visual Studio 2022 using the C++ language.

The modified algorithms (MHB, MNAG) ensure a monotonic decrease of the objective function and prevent oscillations. These modifications incorporate modern line search procedures and restarts. Additionally, we explore the potential for further acceleration of oscillation-free momentum minimizers. Notably, the variation in the friction-related coefficient within the model significantly impacts performance time.