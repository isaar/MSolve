using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Matrices;
using System.Diagnostics;
using System.IO;

namespace ISAAR.MSolve.Solvers.Skyline
{
    //public class ScientificDouble
    //{
    //    public int Exponent { get; set; }
    //    public double Mantissa { get; set; }
    //    public bool IsNegative { get; set; }

    //    public void ReduceAccuracy(int accuracy)
    //    {
    //        if (accuracy < 1) return;

    //        Mantissa = Math.Truncate(Mantissa * Math.Pow(10, accuracy - 1)) * Math.Pow(10, 1 - accuracy);
    //    }
        
    //    public static ScientificDouble GetScientificDouble(double number)
    //    {
    //        ScientificDouble s = new ScientificDouble();
    //        if (number < 0) s.IsNegative = true;
    //        double d = Math.Abs(number);
    //        while (d >= 1)
    //        {
    //            d /= 10;
    //            s.Exponent++;
    //        }
    //        while (d < 1)
    //        {
    //            d *= 10;
    //            s.Exponent--;
    //        }
    //        s.Mantissa = d;
    //        return s;
    //    }

    //    public static double GetDouble(ScientificDouble s)
    //    {
    //        double d = s.Mantissa * Math.Pow(10, s.Exponent);
    //        if (s.IsNegative)
    //            return -1.0 * d;
    //        else
    //            return d;
    //    }
    //}

    public class SolverSkyline : ISolver
    {
        private string stringFormat;
        private int accuracyDigits;
        private readonly Model model;
        private readonly Dictionary<int, ISolverSubdomain> subdomainsDictionary;
        public int AccuracyDigits { get { return accuracyDigits; }
            set 
            { 
                accuracyDigits = value; 
                stringFormat = "";
                if (accuracyDigits < 1) return;

                for (int i = 0; i < accuracyDigits; i++) stringFormat += "#";
                if (accuracyDigits > 1) stringFormat.Insert(1, ".");
                stringFormat += "e+00";
            }
        }

        public SolverSkyline(Model model)
        {
            this.model = model;
            subdomainsDictionary = new Dictionary<int, ISolverSubdomain>(model.SubdomainsDictionary.Count);
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
                subdomainsDictionary.Add(subdomain.ID, new SubdomainSkyline(subdomain));
            this.AccuracyDigits = -1;
        }

        public Dictionary<int, ISolverSubdomain> SubdomainsDictionary
        {
            get { return subdomainsDictionary; }
        }

        #region ISolver Members

        public void Initialize()
        {
            if (model.SubdomainsDictionary.Count != 1) throw new InvalidOperationException("Skyline solver operates on one subdomain only.");
            foreach (ISolverSubdomain subdomain in subdomainsDictionary.Values)
            {
                if (((SkylineMatrix2D<double>)subdomain.Matrix).IsFactorized) continue;

                List<Vector<double>> zems = new List<Vector<double>>();
                List<int> zemColumns = new List<int>();
                SkylineMatrix2D<double> m = (SkylineMatrix2D<double>)subdomain.Matrix;

                //StreamWriter sw = File.CreateText(@"d:\KIx1M.txt");
                //for (int i = 0; i < m.RowIndex.Length; i++)
                //    sw.WriteLine(m.RowIndex[i]);
                //sw.Close();
                //sw = File.CreateText(@"d:\KData1M.txt");
                //for (int i = 0; i < m.Data.Length; i++)
                //    sw.WriteLine(m.Data[i]);
                //sw.Close();
                //sw = File.CreateText(@"d:\f1M.txt");
                //for (int i = 0; i < subdomain.RHS.Length; i++)
                //    sw.WriteLine(subdomain.RHS[i]);
                //sw.Close();
                //return;

                Stopwatch stopWatch = new Stopwatch();
                stopWatch.Start();
                ((SkylineMatrix2D<double>)subdomain.Matrix).Factorize(1e-15, zems, zemColumns);
                stopWatch.Stop();
                //StreamWriter sw = File.CreateText(@"c:\factorization.txt");
                //sw.WriteLine(stopWatch.ElapsedMilliseconds.ToString());
                //sw.Close();
                if (zemColumns.Count > 0) throw new InvalidOperationException("Skyline solver does not operate on singular matrices.");
            }

            //using (var s = File.CreateText("dataFac.txt"))
            //{
            //    var m = (SkylineMatrix2D<double>)subdomainsDictionary[1].Matrix;
            //    for (int i = 0; i < m.Data.Length; i++)
            //        s.WriteLine(String.Format("{0:0.#################E+00}", m.Data[i]));
            //}
        }

        public void LessenAccuracy(double tolerance)
        {
            var valueDictionary = new SortedDictionary<double, List<int>>();
            SkylineMatrix2D<double> k = ((SkylineMatrix2D<double>)subdomainsDictionary[1].Matrix);
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

        private void DestroyAccuracy(ISolverSubdomain subdomain)
        {
            if (AccuracyDigits < 1) return;

            for (int i = 0; i < subdomain.RHS.Length; i++)
            {
                //ScientificDouble s = ScientificDouble.GetScientificDouble(subdomain.RHS[i]);
                //s.ReduceAccuracy(AccuracyDigits);
                //subdomain.RHS[i] = ScientificDouble.GetDouble(s);
                subdomain.RHS[i] = Double.Parse(String.Format("{0:" + stringFormat + "}", subdomain.RHS[i]));
            }
        }

        public void Solve()
        {
            foreach (ISolverSubdomain subdomain in subdomainsDictionary.Values)
            {
                double[] x = ((Vector<double>)subdomain.Solution).Data;
                SkylineMatrix2D<double> k = ((SkylineMatrix2D<double>)subdomain.Matrix);
                Stopwatch stopWatch = new Stopwatch();
                stopWatch.Start();
                k.Solve(subdomain.RHS, x);
                stopWatch.Stop();
                DestroyAccuracy(subdomain);

                x = new double[k.Rows];
                //LessenAccuracy(1e-7);
                k.Solve(subdomain.RHS, x);
                var xVec = new Vector<double>(x);
                var y = xVec.Norm;

                //StreamWriter sw = File.CreateText(@"c:\fbsub.txt");
                //sw.WriteLine(stopWatch.ElapsedMilliseconds.ToString());
                //sw.Close();
            }
        }

        #endregion
    }
}
