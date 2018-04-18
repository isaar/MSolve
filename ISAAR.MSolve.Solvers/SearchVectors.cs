using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System.Threading;
using System.Threading.Tasks;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.Solvers
{
    public static class SearchVectors
    {
        public static void CalculateReorthogonalizedSearchVector(IVector z, IVector p, IList<IVector> ps, IList<IVector> qs)
        {
            int procs = VectorExtensions.AffinityCount;

            //z.CopyTo(p.Data, 0);
            Parallel.ForEach(z.PartitionLimits(procs), limit =>
            {
                //Array.Copy(z.Data, limit.Item2, p.Data, limit.Item2, limit.Item3 - limit.Item2);
                z.CopyFrom(limit.Item2, limit.Item3 - limit.Item2, p, limit.Item2);
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
                            bsTemp1[limit.Item1] += z[i] * qs[c][i];
                            bsTemp2[limit.Item1] += ps[c][i] * qs[c][i];
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
                        for (int j = limit.Item2; j < limit.Item3; j++) p[j] -= ab * ps[c][j];
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

        public static double CalculateReorthogonalizedGradient(IVector p, IVector q, IVector r, IList<IVector> ps, IList<IVector> qs)
        {
            ps.Add(new Vector(p.Length));
            qs.Add(new Vector(q.Length));
            //p.CopyTo(ps[ps.Count - 1].Data, 0);
            //q.CopyTo(qs[qs.Count - 1].Data, 0);
            //return (p * r) / (p * q);

            int procs = VectorExtensions.AffinityCount;
            double[] aTemp = new double[procs];
            double[] bTemp = new double[procs];
            Parallel.ForEach(p.PartitionLimits(procs), limit =>
            {
                //Array.Copy(p.Data, limit.Item2, ps[ps.Count - 1].Data, limit.Item2, limit.Item3 - limit.Item2);
                //Array.Copy(q.Data, limit.Item2, qs[qs.Count - 1].Data, limit.Item2, limit.Item3 - limit.Item2);
                p.CopyFrom(limit.Item2, limit.Item3 - limit.Item2, ps[ps.Count - 1], limit.Item2);
                q.CopyFrom(limit.Item2, limit.Item3 - limit.Item2, qs[ps.Count - 1], limit.Item2);
                for (int j = limit.Item2; j < limit.Item3; j++)
                {
                    aTemp[limit.Item1] += p[j] * r[j];
                    bTemp[limit.Item1] += p[j] * q[j];
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

        public static bool InitializeStartingVectorFromReorthogonalizedSearchVectors(IVector x, IVector b, IList<IVector> ps, IList<IVector> qs)
        {
            if (ps.Count == 0) return false;

            ps.RemoveAt(ps.Count - 1);
            qs.RemoveAt(qs.Count - 1);
            var xx = x;
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
                        aTemp[limit.Item1] += ps[c][i] * qs[c][i];
                        bTemp[limit.Item1] += ps[c][i] * b[i];
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
                    for (int j = limit.Item2; j < limit.Item3; j++)
                        xx[j] += v * ps[c][j];
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
