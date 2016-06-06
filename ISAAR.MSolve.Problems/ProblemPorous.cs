using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.PreProcessor.Providers;
using ISAAR.MSolve.Analyzers;

namespace ISAAR.MSolve.Problems
{
    public class ProblemPorous : IImplicitIntegrationProvider, IStaticProvider, INonLinearProvider
    {
        private bool providersInitialized = false;
        private double scalingCoefficient;
        private Dictionary<int, IMatrix2D<double>> ms, cs, ks;
        private Dictionary<int, Sparse2D<double>> qs;
        private readonly Model model;
        private IDictionary<int, ISolverSubdomain> subdomains;
        private ElementPoreStiffnessProvider stiffnessProvider;
        private ElementPoreDampingProvider dampingProvider;
        private ElementPoreMassProvider massProvider;

        public ProblemPorous(Model model, IDictionary<int, ISolverSubdomain> subdomains)
        {
            this.model = model;
            this.subdomains = subdomains;
        }

        private IDictionary<int, IMatrix2D<double>> Ms
        {
            get
            {
                if (ms == null) BuildMs();
                return ms;
            }
        }

        private IDictionary<int, IMatrix2D<double>> Cs
        {
            get
            {
                if (cs == null) BuildCs();
                return cs;
            }
        }

        private IDictionary<int, IMatrix2D<double>> Ks
        {
            get
            {
                if (ks == null)
                    BuildKs();
                else
                    RebuildKs();
                return ks;
            }
        }

        private void InitializeProvidersAndBuildQs(ImplicitIntegrationCoefficients c)
        {
            //if (providersInitialized) return;

            //providersInitialized = true;
            massProvider = new ElementPoreMassProvider(new ElementStructuralMassProvider(), c.Damping);
            dampingProvider = new ElementPoreDampingProvider(new ElementStructuralDampingProvider(), c.Damping);
            stiffnessProvider = new ElementPoreStiffnessProvider(new ElementStructuralStiffnessProvider(), c.Damping);
            BuildQs();
        }

        private void BuildKs()
        {
            ks = new Dictionary<int, IMatrix2D<double>>(model.SubdomainsDictionary.Count);
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
                ks.Add(subdomain.ID, GlobalMatrixAssemblerSkyline.CalculateGlobalMatrix(subdomain, stiffnessProvider));
        }

        private void RebuildKs()
        {
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                if (subdomain.MaterialsModified)
                    ks[subdomain.ID] = GlobalMatrixAssemblerSkyline.CalculateGlobalMatrix(subdomain, stiffnessProvider);
                subdomain.ResetMaterialsModifiedProperty();
            }
        }

        private void BuildMs()
        {
            ms = new Dictionary<int, IMatrix2D<double>>(model.SubdomainsDictionary.Count);
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
                ms.Add(subdomain.ID, GlobalMatrixAssemblerSkyline.CalculateGlobalMatrix(subdomain, massProvider));
        }

        private void BuildCs()
        {
            cs = new Dictionary<int, IMatrix2D<double>>(model.SubdomainsDictionary.Count);
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
                cs.Add(subdomain.ID, GlobalMatrixAssemblerSkyline.CalculateGlobalMatrix(subdomain, dampingProvider));
        }

        private void BuildQs()
        {
            qs = new Dictionary<int, Sparse2D<double>>(model.SubdomainsDictionary.Count);
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
                qs.Add(subdomain.ID, BuildQFromSubdomain(subdomain));
        }

        private Sparse2D<double> BuildQFromSubdomain(Subdomain subdomain)
        {
            Sparse2D<double> qSubdomain = new Sparse2D<double>(subdomain.TotalDOFs, subdomain.TotalDOFs);
            foreach (Element element in subdomain.ElementsDictionary.Values)
            {
                if (!(element.ElementType is IPorousFiniteElement)) continue;

                IPorousFiniteElement e = (IPorousFiniteElement)element.ElementType;
                IMatrix2D<double> q = e.CouplingMatrix(element);

                int iElementMatrixRow = 0;
                for (int i = 0; i < element.ElementType.DOFEnumerator.GetDOFTypes(element).Count; i++)
                {
                    Node nodeRow = element.Nodes[i];
                    foreach (DOFType dofTypeRow in element.ElementType.DOFEnumerator.GetDOFTypes(element)[i])
                    {
                        if (dofTypeRow != DOFType.Pore) continue;

                        int dofRow = subdomain.NodalDOFsDictionary[nodeRow.ID][dofTypeRow];
                        if (dofRow != -1)
                        {
                            int iElementMatrixColumn = 0;
                            for (int j = 0; j < element.ElementType.DOFEnumerator.GetDOFTypes(element).Count; j++)
                            {
                                Node nodeColumn = element.Nodes[j];
                                foreach (DOFType dofTypeColumn in element.ElementType.DOFEnumerator.GetDOFTypes(element)[j])
                                {
                                    if (dofTypeColumn == DOFType.Pore) continue;

                                    int dofColumn = subdomain.NodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
                                    if (dofColumn != -1)
                                    {
                                        qSubdomain[dofColumn, dofRow] += q[iElementMatrixRow, iElementMatrixColumn];
                                    }
                                    iElementMatrixColumn++;
                                }
                            }
                        }
                        iElementMatrixRow++;
                    }
                }
            }

            return qSubdomain;
        }

        private void ScaleSubdomainSolidVector(ISolverSubdomain subdomain, IVector<double> vector)
        {
            foreach (int nodeID in model.SubdomainsDictionary[subdomain.ID].NodalDOFsDictionary.Keys)
                foreach (DOFType dofType in model.SubdomainsDictionary[subdomain.ID].NodalDOFsDictionary[nodeID].Keys)
                    if (dofType != DOFType.Pore)
                    {
                        int dof = model.SubdomainsDictionary[subdomain.ID].NodalDOFsDictionary[nodeID][dofType];
                        if (dof > -1) vector[dof] *= this.scalingCoefficient;
                    }
        }

        private void CalculateEffectiveMatrixInternal(ISolverSubdomain subdomain, ImplicitIntegrationCoefficients coefficients)
        {
            subdomain.Matrix = this.Ks[subdomain.ID];
            subdomain.Matrix.LinearCombination(
                new double[] 
                {
                    coefficients.Stiffness, coefficients.Mass, coefficients.Damping
                },
                new IMatrix2D<double>[] 
                { 
                    this.Ks[subdomain.ID], this.Ms[subdomain.ID], this.Cs[subdomain.ID] 
                });
        }

        #region IAnalyzerProvider Members
        public void Reset()
        {
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
                foreach (var element in subdomain.ElementsDictionary.Values)
                    element.ElementType.ClearMaterialState();

            cs = null;
            ks = null;
            ms = null;
        }
        #endregion

        #region IImplicitIntegrationProvider Members

        public void CalculateEffectiveMatrix(ISolverSubdomain subdomain, ImplicitIntegrationCoefficients coefficients)
        {
            InitializeProvidersAndBuildQs(coefficients);
            scalingCoefficient = coefficients.Damping;
            CalculateEffectiveMatrixInternal(subdomain, coefficients);
        }

        public void ProcessRHS(ISolverSubdomain subdomain, ImplicitIntegrationCoefficients coefficients)
        {
            scalingCoefficient = coefficients.Damping;
            ScaleSubdomainSolidVector(subdomain, subdomain.RHS);
        }

        public void GetRHSFromHistoryLoad(int timeStep)
        {
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
                for (int i = 0; i < subdomain.Forces.Length; i++)
                    subdomain.Forces[i] = 0;

            model.AssignLoads();
            model.AssignMassAccelerationHistoryLoads(timeStep);
        }

        public void MassMatrixVectorProduct(ISolverSubdomain subdomain, IVector<double> vIn, double[] vOut)
        {
            this.Ms[subdomain.ID].Multiply(vIn, vOut);
        }

        public void DampingMatrixVectorProduct(ISolverSubdomain subdomain, IVector<double> vIn, double[] vOut)
        {
            this.Cs[subdomain.ID].Multiply(vIn, vOut);
            double[] vQ = new double[vOut.Length];
            qs[subdomain.ID].Multiply(vIn, vQ);
            for (int i = 0; i < vOut.Length; i++) vOut[i] -= vQ[i];
        }

        #endregion

        #region IStaticProvider Members

        public void CalculateMatrix(ISolverSubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        //public void CalculateMatrices()
        //{
        //    throw new NotImplementedException();
        //}

        #endregion

        #region INonLinearProvider Members

        public double RHSNorm(double[] rhs)
        {
            //return (new Vector<double>(rhs)).Norm;

            double norm = 0;
            foreach (int nodeID in model.NodalDOFsDictionary.Keys)
                foreach (DOFType dofType in model.NodalDOFsDictionary[nodeID].Keys)
                    if (dofType != DOFType.Pore)
                    {
                        int dof = model.NodalDOFsDictionary[nodeID][dofType];
                        if (dof > -1) norm += rhs[dof] * rhs[dof];
                    }
            return Math.Sqrt(norm);
        }

        public void ProcessInternalRHS(ISolverSubdomain subdomain, double[] rhs, double[] solution)
        {
            //ScaleSubdomainSolidVector(subdomain, new Vector<double>(rhs));
            //return;

            double[] vQ = new double[rhs.Length];
            qs[subdomain.ID].Multiply(new Vector<double>(solution), vQ);
            for (int i = 0; i < rhs.Length; i++) rhs[i] += vQ[i];
            ScaleSubdomainSolidVector(subdomain, new Vector<double>(rhs));
        }

        #endregion
    }
}
