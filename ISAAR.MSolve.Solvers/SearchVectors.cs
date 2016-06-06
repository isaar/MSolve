using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices;
using System.Threading;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Solvers
{
    public static class SearchVectors
    {
        public static void CalculateReorthogonalizedSearchVector(Vector<double> z, Vector<double> p, IList<Vector<double>> ps,
            IList<Vector<double>> qs)
        {
            int procs = VectorExtensions.AffinityCount;

            //z.CopyTo(p.Data, 0);
            Parallel.ForEach(z.PartitionLimits(procs), limit =>
            {
                Array.Copy(z.Data, limit.Item2, p.Data, limit.Item2, limit.Item3 - limit.Item2);
            });
            if (ps.Count > 0)
            {
                //double[] bs = new double[ps.Count];
                double[] bsTemp1 = new double[procs];
                double[] bsTemp2 = new double[procs];
                for (int c = 0; c < ps.Count; c++)
                {
                    Array.Clear(bsTemp1, 0, procs);
                    Array.Clear(bsTemp2, 0, procs);
                    Parallel.ForEach(ps[c].PartitionLimits(procs), limit =>
                    {
                        for (int i = limit.Item2; i < limit.Item3; i++)
                        {
                            bsTemp1[limit.Item1] += z[i] * qs[c].Data[i];
                            bsTemp2[limit.Item1] += ps[c].Data[i] * qs[c].Data[i];
                        }
                    });

                    double a = 0;
                    double b = 0;
                    for (int i = 0; i < procs; i++)
                    {
                        a += bsTemp1[i];
                        b += bsTemp2[i];
                    }
                    double ab = a / b;

                    if (Double.IsNaN(ab)) throw new InvalidOperationException("Gradient of search vector is NaN.");
                    Parallel.ForEach(ps[c].PartitionLimits(procs), limit =>
                    {
                        for (int j = limit.Item2; j < limit.Item3; j++) p.Data[j] -= ab * ps[c].Data[j];
                    });
                }

                //Parallel.For(0, ps.Count, i =>
                //    {
                //        bs[i] = (z * qs[i]) / (ps[i] * qs[i]);
                //    });

                //for (int i = 0; i < ps.Count; i++)
                //{
                //    if (Double.IsNaN(bs[i])) throw new InvalidOperationException("Gradient of search vector is NaN.");
                //    //for (int j = 0; j < p.Length; j++) p.Data[j] -= bs[i] * ps[i].Data[j];
                //    //Parallel.For(0, p.Length, j => { p.Data[j] -= bs[i] * ps[i].Data[j]; });
                //    Parallel.ForEach(ps[i].PartitionLimits(procs), limit =>
                //    {
                //        for (int j = limit.Item2; j < limit.Item3; j++) p.Data[j] -= bs[i] * ps[i].Data[j];
                //    });

                //}
                //for (int i = 0; i < ps.Count; i++)
                //{
                //    double b = (z * qs[i]) / (ps[i] * qs[i]);
                //    if (Double.IsNaN(b)) throw new InvalidOperationException("Gradient of search vector is NaN.");
                //    for (int j = 0; j < p.Length; j++) p.Data[j] -= b * ps[i].Data[j];
                //}
            }
        }

        public static double CalculateReorthogonalizedGradient(Vector<double> p, Vector<double> q, 
            Vector<double> r, IList<Vector<double>> ps, IList<Vector<double>> qs)
        {
            ps.Add(new Vector<double>(p.Length));
            qs.Add(new Vector<double>(q.Length));
            //p.CopyTo(ps[ps.Count - 1].Data, 0);
            //q.CopyTo(qs[qs.Count - 1].Data, 0);
            //return (p * r) / (p * q);

            int procs = VectorExtensions.AffinityCount;
            double[] aTemp = new double[procs];
            double[] bTemp = new double[procs];
            Parallel.ForEach(p.PartitionLimits(procs), limit =>
            {
                Array.Copy(p.Data, limit.Item2, ps[ps.Count - 1].Data, limit.Item2, limit.Item3 - limit.Item2);
                Array.Copy(q.Data, limit.Item2, qs[qs.Count - 1].Data, limit.Item2, limit.Item3 - limit.Item2);
                for (int j = limit.Item2; j < limit.Item3; j++)
                {
                    aTemp[limit.Item1] += p.Data[j] * r.Data[j];
                    bTemp[limit.Item1] += p.Data[j] * q.Data[j];
                }
            });

            double a = 0;
            double b = 0;
            for (int i = 0; i < procs; i++)
            {
                a += aTemp[i];
                b += bTemp[i];
            }
            return a / b;
        }

        public static bool InitializeStartingVectorFromReorthogonalizedSearchVectors(
            Vector<double> x, Vector<double> b, IList<Vector<double>> ps, IList<Vector<double>> qs)
        {
            if (ps.Count == 0) return false;

            ps.RemoveAt(ps.Count - 1);
            qs.RemoveAt(qs.Count - 1);
            Vector<double> xx = (Vector<double>)x;
            xx.Clear();
            double v;

            int procs = VectorExtensions.AffinityCount;
            double[] aTemp = new double[procs];
            double[] bTemp = new double[procs];
            for (int c = 0; c < ps.Count; c++)
            {
                Array.Clear(aTemp, 0, procs);
                Array.Clear(bTemp, 0, procs);
                Parallel.ForEach(ps[c].PartitionLimits(procs), limit =>
                {
                    for (int i = limit.Item2; i < limit.Item3; i++)
                    {
                        aTemp[limit.Item1] += ps[c].Data[i] * qs[c].Data[i];
                        bTemp[limit.Item1] += ps[c].Data[i] * b[i];
                    }
                });

                double a = 0;
                double bT = 0;
                for (int i = 0; i < procs; i++)
                {
                    a += aTemp[i];
                    bT += bTemp[i];
                }
                v = bT / a;

                if (Double.IsNaN(v)) throw new InvalidOperationException("Gradient of starting vector initialization is NaN.");
                Parallel.ForEach(ps[c].PartitionLimits(procs), limit =>
                {
                    for (int j = limit.Item2; j < limit.Item3; j++) xx.Data[j] += v * ps[c].Data[j];
                });
            }
            //for (int i = 0; i < ps.Count; i++)
            //{
            //    v = 1 / (ps[i] * qs[i]) * (ps[i] * (Vector<double>)b);
            //    if (Double.IsNaN(v)) throw new InvalidOperationException("Gradient of starting vector initialization is NaN.");
            //    for (int j = 0; j < xx.Length; j++) xx.Data[j] += v * ps[i].Data[j];
            //}


            return true;
        }
    }
}
