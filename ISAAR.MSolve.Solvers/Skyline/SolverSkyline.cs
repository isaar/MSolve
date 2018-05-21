using System;
using System.Collections.Generic;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Skyline
{
    public class SolverSkyline : ISolver
    {
        private string stringFormat;
        private int accuracyDigits;
        private readonly ILinearSystem linearSystem;

        public SolverSkyline(ILinearSystem linearSystem)
        {
            this.linearSystem = linearSystem;
        }

        #region ISolver Members

        public void Initialize()
        {
            ILinearSystem subdomain = linearSystem;
            if (((SkylineMatrix2D)subdomain.Matrix).IsFactorized) return;

            List<IVector> zems = new List<IVector>();
            List<int> zemColumns = new List<int>();
            SkylineMatrix2D m = (SkylineMatrix2D)subdomain.Matrix;

            ((SkylineMatrix2D)subdomain.Matrix).Factorize(1e-5, zems, zemColumns);
            if (zemColumns.Count > 0) throw new InvalidOperationException("Skyline solver does not operate on singular matrices.");
        }

        public void LessenAccuracy(double tolerance)
        {
            var valueDictionary = new SortedDictionary<double, List<int>>();
            SkylineMatrix2D k = ((SkylineMatrix2D)linearSystem.Matrix);
            for (int i = 0; i < k.Data.Length; i++)
            {
                double difference = Double.MaxValue;
                double targetValue = k.Data[i];
                int index = 0;
                int count = 0;
                foreach (var v in valueDictionary.Keys)
                {
                    if (Math.Abs(k.Data[i] - v) < difference)
                    {
                        difference = Math.Abs(k.Data[i] - v);
                        // Minimize error by taking mean value of these values?
                        targetValue = v;
                        index = count;
                    }
                    count++;
                }
                if (difference > tolerance)
                    valueDictionary.Add(k.Data[i], new List<int>(new[] { i }));
                else
                    valueDictionary[targetValue].Add(i);
            }

            for (int i = 0; i < k.Data.Length; i++)
                foreach (var v in valueDictionary.Keys)
                    if (Math.Abs(k.Data[i] - v) < tolerance)
                    {
                        k.Data[i] = v;
                        break;
                    }
        }

        //private void DestroyAccuracy(ILinearSystem subdomain)
        //{
        //    if (AccuracyDigits < 1) return;

        //    for (int i = 0; i < subdomain.RHS.Length; i++)
        //    {
        //        //ScientificDouble s = ScientificDouble.GetScientificDouble(subdomain.RHS[i]);
        //        //s.ReduceAccuracy(AccuracyDigits);
        //        //subdomain.RHS[i] = ScientificDouble.GetDouble(s);
        //        subdomain.RHS[i] = Double.Parse(String.Format("{0:" + stringFormat + "}", subdomain.RHS[i]));
        //    }
        //}

        public void Solve()
        {
            SkylineMatrix2D k = ((SkylineMatrix2D)linearSystem.Matrix);
            k.Solve(linearSystem.RHS, linearSystem.Solution);
        }

        #endregion
    }
}
