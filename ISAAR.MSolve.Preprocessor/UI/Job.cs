﻿using ISAAR.MSolve.Analyzers;
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
using System.Linq;
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
        /// Defines how loads, boundary conditions and body forces will be incremented in case of material nonlinearity or 
        /// geometric nonlinearity (or both).
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
        /// Defines the linear system solution algorithm.
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

        private readonly Model model;

        /// <summary>
        /// Instantiates an new <see cref="Job"/>. Use its properties to set up the simulation.
        /// </summary>
        /// <param name="model">The mesh, supports, loads, etc. that will be simulated.</param>
        public Job(Model model)
        {
            this.model = model;
        }

        /// <summary>
        /// Defines how loads, boundary conditions and body forces will be incremented in case of material nonlinearity or 
        /// geometric nonlinearity (or both).
        /// </summary>
        public IntegratorOptions Integrator { get; set; } = IntegratorOptions.Linear;

        /// <summary>
        /// Field output requests. It can be left null, in which case, nothing will be output.
        /// </summary>
        public OutputRequests FieldOutputRequests { get; set; } = null;

        /// <summary>
        /// Defines what analysis will be carried out.
        /// </summary>
        public ProcedureOptions Procedure { get; set; } = ProcedureOptions.Static;

        /// <summary>
        /// Defines the linear system solution algorithm.
        /// </summary>
        public SolverOptions Solver { get; set; } = SolverOptions.DirectSkyline;

        public int TimeStep { get; set; } = 0;
        public int TotalTime { get; set; } = 0;

        /// <summary>
        /// Sets up the necessary MSolve objects, checks user input and finally runs the simulation. 
        /// </summary>
        public void Submit()
        {
            // Linear system solver
            VectorExtensions.AssignTotalAffinityCount();
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
                var linearAnalyzer = new LinearAnalyzer(solver, linearSystems);

                // Field output requests 
                //TODO: this should work for all analyzers
                if (FieldOutputRequests != null) linearAnalyzer.LogFactories[0] = FieldOutputRequests.CreateLogFactory(model);

                childAnalyzer = linearAnalyzer;
            }
            else if (Integrator == IntegratorOptions.LoadControl)
            {
                throw new NotImplementedException("The next code doesn't compile since the classes NonLinearSubdomainUpdater"
                    + " and SubdomainGlobalMapping do not exist");
                //INonLinearSubdomainUpdater[] subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
                //ISubdomainGlobalMapping[] subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
                //int increments = 1;
                //int maxIterations = 100;
                //int iterationsForMatrixRebuild = 1;
                //var nonLinearAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystems.Values.ToArray(), 
                //    subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs, maxIterations, 
                //    iterationsForMatrixRebuild);

                //childAnalyzer = nonLinearAnalyzer;
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

            // Run the analysis
            model.ConnectDataStructures();
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }
    }
}
