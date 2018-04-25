using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Solvers.PCG;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.Solvers.PCGSkyline
{
    public class SolverPCG<T> : IIterativeSolver where T : IMatrix2D
    {
        private readonly ILinearSystem linearSystem;
        private readonly SolverPCG solverPCG;
        private readonly ISolverPCGMatrixCalculator matrixCalculator;
        private readonly ISolverPCGInitialization matrixInitialization;
        private readonly ISearchVectorCalculator searchVectorCalculator;
        //private Dictionary<int, Vector<double>> diagonalPreconditioner;
        private IVector x, r;
        private double detf = 0;

        public IVector VectorX { get { return solverPCG.VectorX; } }

        public SolverPCG(ILinearSystem linearSystem, ISearchVectorCalculator searchVectorCalculator, ISolverPCGMatrixCalculator matrixCalculator, ISolverPCGInitialization matrixInitialization)
        {
            this.linearSystem = linearSystem;
            this.matrixCalculator = matrixCalculator;
            this.matrixInitialization = matrixInitialization;
            this.searchVectorCalculator = searchVectorCalculator;
            solverPCG = new SolverPCG(matrixCalculator, searchVectorCalculator);
        }

        public SolverPCG(ILinearSystem linearSystem, ISearchVectorCalculator searchVectorCalculator, ISolverPCGMatrixCalculator matrixCalculator)
        {
            this.linearSystem = linearSystem;
            this.matrixCalculator = matrixCalculator;
            this.matrixInitialization = new SolverPCGMatrixInitialization<T>(this);
            this.searchVectorCalculator = searchVectorCalculator;
            solverPCG = new SolverPCG(matrixCalculator, searchVectorCalculator);
        }

        public SolverPCG(ILinearSystem linearSystem, ISearchVectorCalculator searchVectorCalculator)
        {
            this.linearSystem = linearSystem;
            this.matrixCalculator = new SolverPCGMatrixCalculator<T>(this);
            this.matrixInitialization = new SolverPCGMatrixInitialization<T>(this);
            this.searchVectorCalculator = searchVectorCalculator;
            solverPCG = new SolverPCG(matrixCalculator, searchVectorCalculator);
        }

        public ILinearSystem LinearSystem
        {
            get { return linearSystem; }
        }

        #region ISolver Members

        public void Initialize()
        {
            //if (model.SubdomainsDictionary.Count != 1) throw new InvalidOperationException("Skyline PCG solver operates on one subdomain only.");

            //int vectorLength = matrixCalculator.VectorSize;
            //foreach (ISolverSubdomain subdomain in subdomainsDictionary.Values)
            //{
            //    x = new Vector<double>(vectorLength);
            //    r = new Vector<double>(vectorLength);
            //    SkylineMatrix2D<double> k = ((SkylineMatrix2D<double>)subdomain.Matrix);
            //    diagonalPreconditioner = new Vector<double>(k.Rows);
            //    for (int i = 0; i < subdomain.RHS.Length; i++)
            //    {
            //        detf += subdomain.RHS[i] * subdomain.RHS[i];
            //        r[i] = subdomain.RHS[i];
            //        diagonalPreconditioner[i] = 1 / k.Data[k.RowIndex[i]];
            //    }
            //}
            //detf = Math.Sqrt(detf);

            int vectorLength = matrixCalculator.VectorSize;
            x = new Vector(vectorLength);
            r = new Vector(vectorLength);
            detf = matrixInitialization.InitializeAndGetResidual(new[] { linearSystem }, r, x);
            solverPCG.Initialize(x, r, detf);
        }

        public void Solve()
        {
            Solve(20000, 1e-7);

            linearSystem.Solution = x;
        }

        //public void Solve()
        //{
        //    foreach (ISolverSubdomain subdomain in subdomainsDictionary.Values)
        //    {
        //        double[] x = ((Vector<double>)subdomain.Solution).Data;
        //        SkylineMatrix2D<double> k = ((SkylineMatrix2D<double>)subdomain.Matrix);
        //        k.Solve(subdomain.RHS, x);
        //    }
        //}

        #endregion

        #region IIterativeSolver Members

        public int CurrentIteration
        {
            get { return solverPCG.CurrentIteration; }
        }

        public void Initialize(IVector x0, IVector residual, double detf)
        {
            //throw new InvalidOperationException("SolverPCG<T> solver cannot use this method.");

            int vectorLength = matrixCalculator.VectorSize;
            if (x == null)
                x = new Vector(vectorLength);
            x.CopyFrom(0, vectorLength, x0, 0);
            //Array.Copy(((Vector<double>)x0).Data, x.Data, vectorLength);

            if (r == null)
                r = new Vector(vectorLength);
            r.CopyFrom(0, vectorLength, residual, 0);
            //Array.Copy(((Vector<double>)residual).Data, r.Data, vectorLength);

            double calculatedDetf = 0;
            double temp = 0;
            for (int i = 0; i < r.Length; i++)
            {
                temp = r[i];
                calculatedDetf += temp * temp;
            }
            calculatedDetf = Math.Sqrt(calculatedDetf);

            solverPCG.Initialize(x, r, calculatedDetf);
        }

        public void Solve(int maxIterations, double tolerance)
        {
            solverPCG.Solve(maxIterations, tolerance);
        }

        #endregion

        public void Precondition(IVector vIn, IVector vOut)
        {
            //matrixCalculator.Precondition(vIn, vOut);
            //for (int i = 0; i < vIn.Length; i++) vOut[i] = diagonalPreconditioner[i] * vIn[i];
            for (int i = 0; i < vIn.Length; i++) vOut[i] = 1/this.linearSystem.Matrix[i,i] * vIn[i];
        }

        public void MultiplyWithMatrix(IVector vIn, IVector vOut)
        {
            //matrixCalculator.MultiplyWithMatrix(vIn, vOut);
            //foreach (ISolverSubdomain subdomain in subdomainsDictionary.Values)
            //{
            //    SkylineMatrix2D<double> k = ((SkylineMatrix2D<double>)subdomain.Matrix);
            //    k.Multiply(vIn, ((Vector<double>)vOut).Data);
            //}
            var a = new double[vIn.Length];
            this.linearSystem.Matrix.Multiply(vIn, a);
            for (int i = 0; i < vIn.Length; i++)
            {
                vOut[i] = a[i];
            }
        }
    }
}
