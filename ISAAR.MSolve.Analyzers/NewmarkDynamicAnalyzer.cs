using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using System.Diagnostics;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    public class NewmarkDynamicAnalyzer : IAnalyzer, INonLinearParentAnalyzer
    {
        private readonly double alpha, delta, timeStep, totalTime;
        private double a0, a1, a2, a3, a4, a5, a6, a7;//, a2a0, a3a0, a4a1, a5a1;
        private Dictionary<int, Vector> rhs = new Dictionary<int, Vector>();
        private Dictionary<int, Vector> uu = new Dictionary<int, Vector>();
        private Dictionary<int, Vector> uum = new Dictionary<int, Vector>();
        private Dictionary<int, Vector> uc = new Dictionary<int, Vector>();
        private Dictionary<int, Vector> ucc = new Dictionary<int, Vector>();
        private Dictionary<int, Vector> u = new Dictionary<int, Vector>();
        private Dictionary<int, Vector> v = new Dictionary<int, Vector>();
        private Dictionary<int, Vector> v1 = new Dictionary<int, Vector>();
        private Dictionary<int, Vector> v2 = new Dictionary<int, Vector>();
        //private readonly Dictionary<int, IImplicitIntegrationAnalyzerLog> resultStorages =
        //    new Dictionary<int, IImplicitIntegrationAnalyzerLog>();
        private readonly Dictionary<int, ImplicitIntegrationAnalyzerLog> resultStorages =
            new Dictionary<int, ImplicitIntegrationAnalyzerLog>();
        private readonly IDictionary<int, ILinearSystem> subdomains;
        private readonly IImplicitIntegrationProvider provider;
        private IAnalyzer childAnalyzer;
        private IAnalyzer parentAnalyzer;

        public NewmarkDynamicAnalyzer(IImplicitIntegrationProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains,
            double alpha, double delta, double timeStep, double totalTime)
        {
            this.provider = provider;
            this.childAnalyzer = embeddedAnalyzer;
            this.subdomains = subdomains;
            this.alpha = alpha;
            this.delta = delta;
            this.timeStep = timeStep;
            this.totalTime = totalTime;
            this.childAnalyzer.ParentAnalyzer = this;
        }

        public Dictionary<int, ImplicitIntegrationAnalyzerLog> ResultStorages { get { return resultStorages; } }

        private void InitializeCoefficients()
        {
            if (delta < 0.5) throw new InvalidOperationException("Newmark delta has to be bigger than 0.5.");
            double aLimit = 0.25 * Math.Pow(0.5 + delta, 2);
            if (alpha < aLimit) throw new InvalidOperationException("Newmark alpha has to be bigger than " + aLimit.ToString() + ".");

            a0 = 1 / (alpha * timeStep * timeStep);
            a1 = delta / (alpha * timeStep);
            a2 = 1 / (alpha * timeStep);
            a3 = 1 / (2 * alpha) - 1;
            a4 = delta / alpha - 1;
            a5 = timeStep * 0.5 * (delta / alpha - 2);
            a6 = timeStep * (1 - delta);
            a7 = delta * timeStep;
        }

        public void ResetSolutionVectors()
        {
            foreach (ILinearSystem subdomain in subdomains.Values)
                subdomain.Solution.Clear();
        }

        private void InitializeInternalVectors()
        {
            uu.Clear();
            uum.Clear();
            uc.Clear();
            ucc.Clear();
            u.Clear();
            v.Clear();
            v1.Clear();
            v2.Clear();
            rhs.Clear();

            foreach (ILinearSystem subdomain in subdomains.Values)
            {
                int dofs = subdomain.RHS.Length;
                uu.Add(subdomain.ID, new Vector(dofs));
                uum.Add(subdomain.ID, new Vector(dofs));
                uc.Add(subdomain.ID, new Vector(dofs));
                ucc.Add(subdomain.ID, new Vector(dofs));
                u.Add(subdomain.ID, new Vector(dofs));
                v.Add(subdomain.ID, new Vector(dofs));
                v1.Add(subdomain.ID, new Vector(dofs));
                v2.Add(subdomain.ID, new Vector(dofs));
                rhs.Add(subdomain.ID, new Vector(dofs));

                // Account for initial conditions coming from a previous solution
                subdomain.Solution.CopyTo(v[subdomain.ID].Data, 0);
            }
        }

        private void InitializeMatrices()
        {
            ImplicitIntegrationCoefficients coeffs = new ImplicitIntegrationCoefficients
            {
                Mass = a0,
                Damping = a1,
                Stiffness = 1
            };
            foreach (ILinearSystem subdomain in subdomains.Values)
                provider.CalculateEffectiveMatrix(subdomain, coeffs);

            var m = (SkylineMatrix2D)subdomains[0].Matrix;//TODO: Subdomain matrices should not be retrieved like that.
            var x = new HashSet<double>();
            int nonZeroCount = 0;
            for (int i = 0; i < m.Data.Length; i++)
            {
                nonZeroCount += m.Data[i] != 0 ? 1 : 0;
                if (x.Contains(m.Data[i]) == false)
                    x.Add(m.Data[i]);
            }
            nonZeroCount += 0;
        }

        private void InitializeRHSs()
        {
            ImplicitIntegrationCoefficients coeffs = new ImplicitIntegrationCoefficients
            {
                Mass = a0,
                Damping = a1,
                Stiffness = 1
            };
            foreach (ILinearSystem subdomain in subdomains.Values)
            {
                provider.ProcessRHS(subdomain, coeffs);
                int dofs = subdomain.RHS.Length;
                for (int i = 0; i < dofs; i++) rhs[subdomain.ID][i] = subdomain.RHS[i];
            }
        }

        private void UpdateResultStorages(DateTime start, DateTime end)
        {
            if (childAnalyzer == null) return;
            foreach (ILinearSystem subdomain in subdomains.Values)
                if (resultStorages.ContainsKey(subdomain.ID))
                    if (resultStorages[subdomain.ID] != null)
                        foreach (var l in childAnalyzer.Logs[subdomain.ID])
                            resultStorages[subdomain.ID].StoreResults(start, end, l);
        }

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog[]> Logs { get { return null; } }
        public IAnalyzer ParentAnalyzer
        {
            get { return parentAnalyzer; }
            set { parentAnalyzer = value; }
        }

        public IAnalyzer ChildAnalyzer
        {
            get { return childAnalyzer; }
            set { childAnalyzer = value; }
        }

        public void Initialize()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");
            //InitializeCoefficients();
            InitializeInternalVectors();
            //InitializeMatrices();
            InitializeRHSs();
            childAnalyzer.Initialize();
        }

        private void CalculateRHSImplicit(ILinearSystem subdomain, double[] rhsResult, bool addRHS)
        {
            int id = subdomain.ID;
            for (int i = 0; i < subdomain.RHS.Length; i++)
            {
                //uu[id].Data[i] = v[id].Data[i] + a2a0 * v1[id].Data[i] + a3a0 * v2[id].Data[i];
                uu[id].Data[i] = a0 * v[id].Data[i] + a2 * v1[id].Data[i] + a3 * v2[id].Data[i];
                uc[id].Data[i] = a1 * v[id].Data[i] + a4 * v1[id].Data[i] + a5 * v2[id].Data[i];
            }
            //provider.Ms[id].Multiply(uu[id], uum[id].Data);
            provider.MassMatrixVectorProduct(subdomain, uu[id], uum[id].Data);
            provider.DampingMatrixVectorProduct(subdomain, uc[id], ucc[id].Data);
            if (addRHS)
                for (int i = 0; i < subdomain.RHS.Length; i++)
                    rhsResult[i] = rhs[id].Data[i] + uum[id].Data[i] + ucc[id].Data[i];
            else
                for (int i = 0; i < subdomain.RHS.Length; i++)
                    rhsResult[i] = uum[id].Data[i] + ucc[id].Data[i];
            #region Old Fortran code
            //If (ptAnal.tDynamic.bHasDamping) Then
            //    atDynUU(iMesh).afValues = atDynV(iMesh).afValues + fA2A0 * atDynV1(iMesh).afValues + fA3A0 * atDynV2(iMesh).afValues
            //    atDynUC(iMesh).afValues = atDynV(iMesh).afValues + fA4A1 * atDynV1(iMesh).afValues + fA5A1 * atDynV2(iMesh).afValues
            //Else
            //    atDynUU(iMesh).afValues = atDynV(iMesh).afValues + fA2A0 * atDynV1(iMesh).afValues + fA3A0 * atDynV2(iMesh).afValues
            //End If
            //If (Present(ProcessVel)) Then
            //    Call ProcessVel(iMesh, Size(atDynUU(iMesh).afValues), atDynUU(iMesh).afValues, atDynUUM(iMesh).afValues)
            //Else
            //    Call MatrixVectorMultiplication(iMesh, Size(atRHS(iMesh).afValues), 2, atDynUU(iMesh).afValues, atDynUUM(iMesh).afValues)
            //End If
            //If (ptAnal.tDynamic.bHasDamping) Then
            //    If (Present(ProcessAcc)) Then
            //        Call ProcessAcc(iMesh, Size(atDynUU(iMesh).afValues), atDynUC(iMesh).afValues, atDynUCC(iMesh).afValues)
            //    Else
            //        Call MatrixVectorMultiplication(iMesh, Size(atRHS(iMesh).afValues), 3, atDynUC(iMesh).afValues, atDynUCC(iMesh).afValues)
            //    End If
            //    atRHS(iMesh).afValues = atRHS(iMesh).afValues + atDynUUM(iMesh).afValues + atDynUCC(iMesh).afValues
            //Else
            //    atRHS(iMesh).afValues = atRHS(iMesh).afValues + atDynUUM(iMesh).afValues
            //End If
            #endregion
        }

        private void CalculateRHSImplicit()
        {
            foreach (ILinearSystem subdomain in subdomains.Values)
            {
                //int id = subdomain.ID;
                //for (int i = 0; i < subdomain.RHS.Length; i++)
                //    subdomain.RHS[i] = rhs[id].Data[i];

                CalculateRHSImplicit(subdomain, ((Vector)subdomain.RHS).Data, true);

                //int id = subdomain.ID;
                //for (int i = 0; i < subdomain.RHS.Length; i++)
                //{
                //    //uu[id].Data[i] = v[id].Data[i] + a2a0 * v1[id].Data[i] + a3a0 * v2[id].Data[i];
                //    uu[id].Data[i] = a0 * v[id].Data[i] + a2 * v1[id].Data[i] + a3 * v2[id].Data[i];
                //    uc[id].Data[i] = a1 * v[id].Data[i] + a4 * v1[id].Data[i] + a5 * v2[id].Data[i];
                //}
                ////provider.Ms[id].Multiply(uu[id], uum[id].Data);
                //provider.MassMatrixVectorProduct(subdomain, uu[id], uum[id].Data);
                //provider.DampingMatrixVectorProduct(subdomain, uc[id], ucc[id].Data);
                //for (int i = 0; i < subdomain.RHS.Length; i++)
                //    subdomain.RHS[i] = rhs[id].Data[i] + uum[id].Data[i] + ucc[id].Data[i];
            }
            #region Old Fortran code
            //If (ptAnal.tDynamic.bHasDamping) Then
            //    atDynUU(iMesh).afValues = atDynV(iMesh).afValues + fA2A0 * atDynV1(iMesh).afValues + fA3A0 * atDynV2(iMesh).afValues
            //    atDynUC(iMesh).afValues = atDynV(iMesh).afValues + fA4A1 * atDynV1(iMesh).afValues + fA5A1 * atDynV2(iMesh).afValues
            //Else
            //    atDynUU(iMesh).afValues = atDynV(iMesh).afValues + fA2A0 * atDynV1(iMesh).afValues + fA3A0 * atDynV2(iMesh).afValues
            //End If
            //If (Present(ProcessVel)) Then
            //    Call ProcessVel(iMesh, Size(atDynUU(iMesh).afValues), atDynUU(iMesh).afValues, atDynUUM(iMesh).afValues)
            //Else
            //    Call MatrixVectorMultiplication(iMesh, Size(atRHS(iMesh).afValues), 2, atDynUU(iMesh).afValues, atDynUUM(iMesh).afValues)
            //End If
            //If (ptAnal.tDynamic.bHasDamping) Then
            //    If (Present(ProcessAcc)) Then
            //        Call ProcessAcc(iMesh, Size(atDynUU(iMesh).afValues), atDynUC(iMesh).afValues, atDynUCC(iMesh).afValues)
            //    Else
            //        Call MatrixVectorMultiplication(iMesh, Size(atRHS(iMesh).afValues), 3, atDynUC(iMesh).afValues, atDynUCC(iMesh).afValues)
            //    End If
            //    atRHS(iMesh).afValues = atRHS(iMesh).afValues + atDynUUM(iMesh).afValues + atDynUCC(iMesh).afValues
            //Else
            //    atRHS(iMesh).afValues = atRHS(iMesh).afValues + atDynUUM(iMesh).afValues
            //End If
            #endregion
        }

        public void Solve()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");

            for (int i = 0; i < (int)(totalTime / timeStep); i++)
            {
                Debug.WriteLine("Newmark step: {0}", i);
                provider.GetRHSFromHistoryLoad(i);
                InitializeRHSs();
                // ProcessRHS
                CalculateRHSImplicit();
                DateTime start = DateTime.Now;
                childAnalyzer.Solve();
                DateTime end = DateTime.Now;
                UpdateVelocityAndAcceleration(i);
                UpdateResultStorages(start, end);
            }

            //childAnalyzer.Solve();
        }

        private void UpdateVelocityAndAcceleration(int timeStep)
        {
            var externalVelocities = provider.GetVelocitiesOfTimeStep(timeStep);
            var externalAccelerations = provider.GetAccelerationsOfTimeStep(timeStep);
            foreach (ILinearSystem subdomain in subdomains.Values)
            {
                int id = subdomain.ID;
                v[id].CopyTo(u[id].Data, 0);
                subdomain.Solution.CopyTo(v[id].Data, 0);
                for (int j = 0; j < subdomain.RHS.Length; j++)
                {
                    double vv = v2[id].Data[j] + externalAccelerations[id][j];
                    v2[id].Data[j] = a0 * (v[id].Data[j] - u[id].Data[j]) - a2 * v1[id].Data[j] - a3 * vv;
                    v1[id].Data[j] += externalVelocities[id][j] + a6 * vv + a7 * v2[id].Data[j];
                }
            }
        }

        public void BuildMatrices()
        {
            InitializeCoefficients();
            InitializeMatrices();
        }

        #endregion

        #region INonLinearParentAnalyzer Members

        public double[] GetOtherRHSComponents(ILinearSystem subdomain, IVector currentSolution)
        {
            //// u[id]: old solution
            //// v[id]: current solution
            //// vv: old acceleration
            //// v2: current acceleration
            //// v1: current velocity
            ////double vv = v2[id].Data[j];
            ////v2[id].Data[j] = a0 * (v[id].Data[j] - u[id].Data[j]) - a2 * v1[id].Data[j] - a3 * vv;
            ////v1[id].Data[j] += a6 * vv + a7 * v2[id].Data[j];

            //int id = subdomain.ID;
            //Vector<double> currentAcceleration = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> currentVelocity = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> uu = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> uc = new Vector<double>(subdomain.Solution.Length);
            //for (int j = 0; j < subdomain.RHS.Length; j++)
            //{
            //    currentAcceleration.Data[j] = a0 * (currentSolution[j] - v[id].Data[j]) - a2 * v1[id].Data[j] - a3 * v2[id].Data[j];
            //    currentVelocity.Data[j] = v1[id].Data[j] + a6 * v2[id].Data[j] + a7 * currentAcceleration.Data[j];
            //    uu.Data[j] = a0 * currentSolution[j] + a2 * currentVelocity.Data[j] + a3 * currentAcceleration.Data[j];
            //    uc.Data[j] = a1 * currentSolution[j] + a4 * currentVelocity.Data[j] + a5 * currentAcceleration.Data[j];
            //}

            //Vector<double> tempResult = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> result = new Vector<double>(subdomain.Solution.Length);
            //provider.MassMatrixVectorProduct(subdomain, uu, tempResult.Data);
            //result.Add(tempResult);

            //provider.DampingMatrixVectorProduct(subdomain, uc, tempResult.Data);
            //result.Add(tempResult);

            //return result.Data;

            Vector result = new Vector(subdomain.Solution.Length);
            Vector temp = new Vector(subdomain.Solution.Length);
            Vector tempResult = new Vector(subdomain.Solution.Length);
            //Vector<double> uu = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> uc = new Vector<double>(subdomain.Solution.Length);
            //int id = subdomain.ID;
            //for (int j = 0; j < subdomain.RHS.Length; j++)
            //{
            //    uu.Data[j] = -a0 * (v[id].Data[j] - currentSolution[j]) - a2 * v1[id].Data[j] - a3 * v2[id].Data[j];
            //    uc.Data[j] = -a1 * (v[id].Data[j] - currentSolution[j]) - a4 * v1[id].Data[j] - a5 * v2[id].Data[j];
            //}
            //provider.MassMatrixVectorProduct(subdomain, uu, tempResult.Data);
            //result.Add(tempResult);
            //provider.DampingMatrixVectorProduct(subdomain, uc, tempResult.Data);
            //result.Add(tempResult);

            ////CalculateRHSImplicit(subdomain, result.Data, false);
            ////result.Scale(-1d);
            currentSolution.CopyTo(temp.Data, 0);
            temp.Scale(a0);
            provider.MassMatrixVectorProduct(subdomain, temp, tempResult.Data);
            result.Add(tempResult);

            currentSolution.CopyTo(temp.Data, 0);
            temp.Scale(a1);
            provider.DampingMatrixVectorProduct(subdomain, temp, tempResult.Data);
            result.Add(tempResult);

            return result.Data;
        }

        #endregion
    }
}
