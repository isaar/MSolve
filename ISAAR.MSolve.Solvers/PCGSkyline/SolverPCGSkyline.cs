using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Solvers.PCG;

namespace ISAAR.MSolve.Solvers.PCGSkyline
{
    public class SolverPCG<T> : IIterativeSolver where T : IMatrix2D<double>
    {
        private readonly Model model;
        private readonly Dictionary<int, ISolverSubdomain> subdomainsDictionary;
        private readonly SolverPCG solverPCG;
        private readonly ISolverPCGMatrixCalculator matrixCalculator;
        private readonly ISolverPCGInitialization matrixInitialization;
        private readonly ISearchVectorCalculator searchVectorCalculator;
        //private Dictionary<int, Vector<double>> diagonalPreconditioner;
        private Vector<double> x, r;
        private double detf = 0;

        public Vector<double> VectorX { get { return solverPCG.VectorX; } }

        public SolverPCG(Model model, ISearchVectorCalculator searchVectorCalculator, ISolverPCGMatrixCalculator matrixCalculator, ISolverPCGInitialization matrixInitialization)
        {
            this.model = model;
            this.matrixCalculator = matrixCalculator;
            this.matrixInitialization = matrixInitialization;
            this.searchVectorCalculator = searchVectorCalculator;
            solverPCG = new SolverPCG(matrixCalculator, searchVectorCalculator);
            subdomainsDictionary = new Dictionary<int, ISolverSubdomain>(model.SubdomainsDictionary.Count);
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
                subdomainsDictionary.Add(subdomain.ID, new SubdomainSkyline(subdomain));
        }

        public SolverPCG(Model model, ISearchVectorCalculator searchVectorCalculator, ISolverPCGMatrixCalculator matrixCalculator)
        {
            this.model = model;
            this.matrixCalculator = matrixCalculator;
            this.matrixInitialization = new SolverPCGMatrixInitialization<T>(this);
            this.searchVectorCalculator = searchVectorCalculator;
            solverPCG = new SolverPCG(matrixCalculator, searchVectorCalculator);
            subdomainsDictionary = new Dictionary<int, ISolverSubdomain>(model.SubdomainsDictionary.Count);
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
                subdomainsDictionary.Add(subdomain.ID, new SubdomainSkyline(subdomain));
        }

        public SolverPCG(Model model, ISearchVectorCalculator searchVectorCalculator)
        {
            this.model = model;
            this.matrixCalculator = new SolverPCGMatrixCalculator<T>(this);
            this.matrixInitialization = new SolverPCGMatrixInitialization<T>(this);
            this.searchVectorCalculator = searchVectorCalculator;
            solverPCG = new SolverPCG(matrixCalculator, searchVectorCalculator);
            subdomainsDictionary = new Dictionary<int, ISolverSubdomain>(model.SubdomainsDictionary.Count);
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
                subdomainsDictionary.Add(subdomain.ID, new SubdomainSkyline(subdomain));
        }

        private void Initialize(Model model, ISearchVectorCalculator searchVectorCalculator, ISolverPCGMatrixCalculator matrixCalculator)
        {
        }

        public Dictionary<int, ISolverSubdomain> SubdomainsDictionary
        {
            get { return subdomainsDictionary; }
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
            x = new Vector<double>(vectorLength);
            r = new Vector<double>(vectorLength);
            detf = matrixInitialization.InitializeAndGetResidual(subdomainsDictionary.Select(s => s.Value).ToArray<ISolverSubdomain>(), r.Data, x.Data);
            solverPCG.Initialize(x, r, detf);
        }

        public void Solve()
        {
            Solve(20000, 1e-7);

            foreach (ISolverSubdomain subdomain in subdomainsDictionary.Values)
                subdomain.Solution = x;
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

        public void Initialize(IVector<double> x0, IVector<double> residual, double detf)
        {
            //throw new InvalidOperationException("SolverPCG<T> solver cannot use this method.");

            int vectorLength = matrixCalculator.VectorSize;
            if (x == null)
                x = new Vector<double>(vectorLength);
            Array.Copy(((Vector<double>)x0).Data, x.Data, vectorLength);

            if (r == null)
                r = new Vector<double>(vectorLength);
            Array.Copy(((Vector<double>)residual).Data, r.Data, vectorLength);

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

        public void Precondition(IVector<double> vIn, IVector<double> vOut)
        {
            matrixCalculator.Precondition(vIn, vOut);
            //for (int i = 0; i < vIn.Length; i++) vOut[i] = diagonalPreconditioner[i] * vIn[i];
        }

        public void MultiplyWithMatrix(IVector<double> vIn, IVector<double> vOut)
        {
            matrixCalculator.MultiplyWithMatrix(vIn, vOut);
            //foreach (ISolverSubdomain subdomain in subdomainsDictionary.Values)
            //{
            //    SkylineMatrix2D<double> k = ((SkylineMatrix2D<double>)subdomain.Matrix);
            //    k.Multiply(vIn, ((Vector<double>)vOut).Data);
            //}
        }
    }
}
