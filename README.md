# Experiments On Numerical Optimization

folder OSC_FREE_MOMENTUMS:
## Oscillation-free Momentum Optimization Algorithms
## არაოსცილირებადი ინერციიანი ოპტიმიზირების ალგორითმები

Folder “experiments-on-numerical-optimization” consists of unconstrained minimization algorithms, and test functions to benchmark them. 
Here, four solvers are considered. They are built in Visual Studio environment using C++.
 
Two of them are well-known as best choice: lcg (W. W. Hager, H. Zhang, The limited memory conjugate gradient method. SIAM J. Optim. 23 (2013), no. 4, 2150–2168) and L-BFGS (Dong C. Liu, Jorge Nocedal, On the limited memory BFGS method for large scale optimization. Math. Programming 45 (1989), no. 3 (Ser. B), 503–528).
 
Project CG_DESCENT-C-6.8 is lcg, adapted to Visual Studio. Code is taken from [http://people.clas.ufl.edu/hager/files/CG_DESCENT-C-6.8.tar_.gz](http://people.clas.ufl.edu/hager/files/CG_DESCENT-C-6.8.tar_.gz).
 
Projects L_BFGS_22 and L-BFGS_LB contain our C++ implementation of the L-BFGS algorithm. The Line Search algorithm is mainly the same in all projects.
 
Two other solvers are Oscillation-free Momentum Optimization Algorithms: modified Polyak algorithm (Shortly MHB) and modified Nesterov algorithm (MNAG). Both in C++.
  
In the MHB, when the value of the objective function starts to increase during the standard HB, the MHB is restarted, guaranteeing reducing the value of the function (at least in the first iteration). At each restart, one of the best line search algorithms is used to increase the initial momentum and to determine the learning rate coefficient. As a result of these changes, MHB became competitive with lcg and L-BFGS, at least for the collection of test functions that were considered in numerical experiments in [Koba Gelashvili, Irina Khutsishvili, Luka Gorgadze, and Lela Alkhazishvili. Speeding up the convergence of the polyak’s heavy ball algorithm. Transactions of A. Razmadze Mathematical Institute, 172 no. 2:176–188, 2018].
 
The generalization of this experiment means that the same can be done for other momentum algorithms. Particularly, the same changes that have been done with HB were added to Nesterov accelerated gradient method, i.e. the Line Search and restarts. As a result, a fast oscillation-free momentum algorithm was obtained.
Projects MHB_1 and MHB_LB contain implementations of MHB, but the newest version of this algorithm is in the project OSC_FREE_MOMENTUMS. 
Projects MHB_LB and L-BFGS_LB  contain add-in  for solvers of unconstrained minimization to eliminate lower bounds of variables by transformation. The add-in is implemented in C++ and consists of three components: small collection of transformations of variables; small collection of minimization test-problems with constraints only on lower bounds of variables; stopping conditions of minimization algorithms and program-driver.

In the numerical experiments, mainly CG-DESCENT-C-6.8 and l-bfgs were used. Results show that empowering unconstrained minimization algorithms by such add-in is a promising direction, especially, for symmetric matrix games. Fully on this subject, see: Koba Gelashvili. Add-in for solvers of unconstrained minimization to eliminate lower bounds of variables by transformation. Transactions of A. Razmadze Mathematical Institute, 173:39–46, 2019.

The project OSC_FREE_MOMENTUMS experiments with the creation of modified, faster momentum minimization algorithms based on existing ones. The modified algorithms monotonically decrease the objective function and does not allow it to oscillate. Such an approach has already been implemented for the classical momentum algorithm, i.e. the Heavy Ball algorithm of the Polyak, and now it is successfully applied to Nesterov accelerated gradient method. Modified Oscillation-free Momentum Optimization Algorithms are equipped with modern line search procedures and restarts.

To determine the efficiency of the new algorithm, numerical experiments were conducted on standard optimization test functions and on single-layer neural networks for several popular datasets. Comparisons were made with the best unconstrained minimization algorithms – lcg, lbgfs and ADAM. Recent experiments have shown that the modified Nesterov accelerated gradient method is even more promising than Modified Heavy Ball, especially in neural networks. It is experimentally proven that the modified algorithms have a specific resource for further increasing the speed. In particular, the value range of one parameter of basic algorithms are wider, which creates new possibilities. 

