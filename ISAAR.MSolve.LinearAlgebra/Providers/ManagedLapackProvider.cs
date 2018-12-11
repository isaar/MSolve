using System;
using System.Collections.Generic;
using System.Text;
using DotNumerics.LinearAlgebra.CSLapack;
using ISAAR.MSolve.LinearAlgebra.Providers.Implementations;

namespace ISAAR.MSolve.LinearAlgebra.Providers
{
    public class ManagedLapackProvider : ILapackProvider
    {
        private static readonly DGELQF dgelqf = new DGELQF();
        private static readonly DGEQRF dgeqrf = new DGEQRF();
        private static readonly DORGLQ dorglq = new DORGLQ();
        private static readonly DORGQR dorgqr = new DORGQR();
        private static readonly DORMLQ dormlq = new DORMLQ();
        private static readonly DORMQR dormqr = new DORMQR();
        private static readonly DGETRF dgetrf = new DGETRF();
        private static readonly DGETRI dgetri = new DGETRI();
        private static readonly DGETRS dgetrs = new DGETRS();

        public void Dgelqf(int m, int n, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, 
            double[] work, int offsetWork, int lWork, ref int info)
            => dgelqf.Run(m, n, ref a, offsetA, ldA, ref tau, offsetTau, ref work, offsetWork, lWork, ref info);

        public void Dgeqrf(int m, int n, double[] a, int offsetA, int ldA, double[] tau, int offsetTau,
            double[] work, int offsetWork, int lWork, ref int info)
            => dgeqrf.Run(m, n, ref a, offsetA, ldA, ref tau, offsetTau, ref work, offsetWork, lWork, ref info);

        public void Dgetrf(int m, int n, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv, ref int info)
            => dgetrf.Run(m, n, ref a, offsetA, ldA, ref ipiv, offsetIpiv, ref info);

        public void Dgetri(int n, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv, 
            double[] work, int offsetWork, int lWork, ref int info)
            => dgetri.Run(n, ref a, offsetA, ldA, ipiv, offsetIpiv, ref work, offsetWork, lWork, ref info);

        public void Dgetrs(string transA, int n, int nRhs, double[] a, int offsetA, int ldA, int[] ipiv, int offsetIpiv,
            double[] b, int offsetB, int ldB, ref int info)
            => dgetrs.Run(transA, n, nRhs, a, offsetA, ldA, ipiv, offsetIpiv, ref b, offsetB, ldB, ref info);

        public void Dorglq(int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, int offsetTau, 
            double[] work, int offsetWork, int lWork, ref int info)
            => dorglq.Run(m, n, k, ref a, offsetA, ldA, tau, offsetTau, ref work, offsetWork, lWork, ref info);

        public void Dorgqr(int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, int offsetTau,
            double[] work, int offsetWork, int lWork, ref int info)
            => dorgqr.Run(m, n, k, ref a, offsetA, ldA, tau, offsetTau, ref work, offsetWork, lWork, ref info);

        public void Dormlq(string side, string transQ, int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau, 
            int offsetTau, double[] c, int offsetC, int ldC, double[] work, int offsetWork, int lWork, ref int info)
            => dormlq.Run(side, transQ, m, n, k, ref a, offsetA, ldA, tau, offsetTau, ref c, offsetC, ldC,
                ref work, offsetWork, lWork, ref info);

        public void Dormqr(string side, string transQ, int m, int n, int k, double[] a, int offsetA, int ldA, double[] tau,
            int offsetTau, double[] c, int offsetC, int ldC, double[] work, int offsetWork, int lWork, ref int info)
            => dormqr.Run(side, transQ, m, n, k, ref a, offsetA, ldA, tau, offsetTau, ref c, offsetC, ldC,
                ref work, offsetWork, lWork, ref info);

        public void Dpotrf(string uplo, int n, double[] a, int offsetA, int ldA, ref int info)
        {
            throw new NotImplementedException();
        }

        public void Dpotri(string uplo, int n, double[] a, int offsetA, int ldA, ref int info)
        {
            throw new NotImplementedException();
        }

        public void Dpotrs(string uplo, int n, int nRhs, double[] a, int offsetA, int ldA, 
            double[] b, int offsetB, int ldB, ref int info)
        {
            throw new NotImplementedException();
        }

        public void Dpptrf(string uplo, int n, double[] a, int offsetA, ref int info)
        {
            throw new NotImplementedException();
        }

        public void Dpptri(string uplo, int n, double[] a, int offsetA, ref int info)
        {
            throw new NotImplementedException();
        }

        public void Dpptrs(string uplo, int n, int nRhs, double[] a, int offsetA, double[] b, int offsetB, int ldB, ref int info)
        {
            throw new NotImplementedException();
        }
    }
}
