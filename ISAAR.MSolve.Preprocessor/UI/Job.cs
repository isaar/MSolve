using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.PCG;
using ISAAR.MSolve.Solvers.PCGSkyline;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Text;

//TODO: figure out how to set up different problem types. Only ProblemStructual is supported now.
//TODO: add non linear and dynamic analyzers
//TODO: provide the option to explicitly set one of the modules, while choosing the presets for the rest.
//TODO: output requests should be appliable to all (child) analyzers
//TODO: time steps and total time should be accessed from the dynamic loads. Ideally the dynamic analyzer itself should access 
//      them from the dynamic loads.
namespace ISAAR.MSolve.Preprocessor.UI
{
    /// <summary>
    /// Utility class for setting up common analyses. Eventually the analyzers should be user-friendly enough to not need this 
    /// class, but at the time this is written, it is useful to abstract the interconnection between analyzers and other classes.
    /// This is bad code. It doesn't follow OOP and violates OCP. Still it is convenient for new users, since it doesn't force 
    /// them to know how every module interacts with the rest.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Job
    {
        /// <summary>
        /// Defines how loads, boundary conditions and body forces will be incremented in case of nonlinear material or 
        /// nonlinear geometry (or both).
        /// </summary>
        public enum IntegratorOptions
        {
            /// <summary>
            /// Loading will be applied in a single step, without nested Newton-Raphson iterations. Appropriate for pure linear
            /// analysis (elastic materials and small deformations).
            /// </summary>
            Linear,

            /// <summary>
            /// Loading will be increased in equal steps. Nested Newton-Rapshon iterations.
            /// </summary>
            LoadControl,

            /// <summary>
            /// Loading will be increased, such that the displacement increments of a sentinel freedom degree are constant.
            /// Nested Newton-Rapshon iterations.
            /// </summary>
            DisplacementControl,

            /// <summary>
            /// TODO: someone else should explain this
            /// </summary>
            ArcLength
        }

        /// <summary>
        /// Defines what analysis will be carried out.
        /// </summary>
        public enum ProcedureOptions
        {
            /// <summary>
            /// 
            /// </summary>
            Static,

            /// <summary>
            /// 
            /// </summary>
            DynamicExplicit,

            /// <summary>
            /// 
            /// </summary>
            DynamicImplicit,

            /// <summary>
            /// 
            /// </summary>
            Modal,

            /// <summary>
            /// 
            /// </summary>
            Buckling,

            /// <summary>
            /// 
            /// </summary>
            MonteCarlo,

            /// <summary>
            /// 
            /// </summary>
            Optimization
        }

        /// <summary>
        /// Defines 
        /// </summary>
        public enum SolverOptions
        {
            /// <summary>
            /// Direct solver using Cholesky factorization on matrices in skyline storage format. It works if there is only one 
            /// subdomain. Higher memory requirements then iterative solvers, but usually faster (especially for 1D and 2D 
            /// problems).
            /// </summary>
            DirectSkyline,

            /// <summary>
            /// Iterative solver using the Preconditioned Conjugate Gradient algorithm on matrices in CSR storage format. It 
            /// works if there is only one subdomain. Use this if the memory requirements for <see cref="DirectSkyline"/> cannot 
            /// be met.
            /// </summary>
            IterativePcg,

            /// <summary>
            /// Domain decomposition solver based on the Finite Element Tearing and Interconnecting (FETI) method.
            /// </summary>
            SubdomainsDual,

            /// <summary>
            /// Domain decomposition solver based on the Primal Substructuring Method.
            /// </summary>
            SubdomainsPrimal            
        }

        public IntegratorOptions Integrator { get; set; } = IntegratorOptions.Linear;
        public ProcedureOptions Procedure { get; set; } = ProcedureOptions.Static;
        public SolverOptions Solver { get; set; } = SolverOptions.DirectSkyline;

        public int TimeStep { get; set; } = 0;
        public int TotalTime { get; set; } = 0;

        public void Submit(Model model, OutputRequests output)
        {
            // Linear system solver
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[0] = new SkylineLinearSystem(0, model.Subdomains[0].Forces);
            ISolver solver;
            if (Solver == SolverOptions.DirectSkyline)
            {
                solver = new SolverSkyline(linearSystems[0]);
            }
            else if (Solver == SolverOptions.IterativePcg)
            {
                solver = new SolverPCG<SkylineMatrix2D>(linearSystems[0], new SolverPCGSimpleSearchVectorCalculator());
            }
            else
            {
                throw new NotImplementedException();
            }

            // Provider of the problem
            ProblemStructural provider = new ProblemStructural(model, linearSystems); //TODO: extend this

            // (Non) linear analyzer
            IAnalyzer childAnalyzer;
            if (Integrator == IntegratorOptions.Linear)
            {
                // Field output requests
                var linearAnalyzer = new LinearAnalyzer(solver, linearSystems);
                linearAnalyzer.LogFactories[0] = output.CreateLogFactory(model);
                childAnalyzer = linearAnalyzer;
            }
            else if (Integrator == IntegratorOptions.LoadControl)
            {
                throw new NotImplementedException();
                //childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystems, )
            }
            else
            {
                throw new NotImplementedException();
            }

            // Parent analyzer
            IAnalyzer parentAnalyzer;
            if (Procedure == ProcedureOptions.Static)
            {
                parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
            }
            else if (Procedure == ProcedureOptions.DynamicImplicit)
            {
                parentAnalyzer = new NewmarkDynamicAnalyzer(provider, childAnalyzer, linearSystems, 0.6, 1, TimeStep, TotalTime);
            }
            else
            {
                throw new NotImplementedException();
            }

            // Field output requests

            // Run the analysis
            model.ConnectDataStructures();
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }
    }
}
