using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Problems
{
    public class ProblemPorous : IImplicitIntegrationProvider, IStaticProvider, INonLinearProvider
    {
        private bool providersInitialized = false;
        private double scalingCoefficient;
        private Dictionary<int, IMatrix> ms, cs, ks;
        private Dictionary<int, CsrMatrix> qs;
        private readonly Model model;
        private readonly ISolver solver;
        private readonly IReadOnlyDictionary<int, ILinearSystem> linearSystems;
        private ElementPoreStiffnessProvider stiffnessProvider;
        private ElementPoreDampingProvider dampingProvider;
        private ElementPoreMassProvider massProvider;

        public ProblemPorous(Model model, ISolver solver)
        {
            this.model = model;
            this.solver = solver;
            this.linearSystems = solver.LinearSystems;
        }

        public IDirichletEquivalentLoadsAssembler DirichletLoadsAssembler => throw new NotImplementedException();

        private IDictionary<int, IMatrix> Ms
        {
            get
            {
                if (ms == null) BuildMs();
                return ms;
            }
        }

        private IDictionary<int, IMatrix> Cs
        {
            get
            {
                if (cs == null) BuildCs();
                return cs;
            }
        }

        private IDictionary<int, IMatrix> Ks
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

        private void BuildKs() => ks = solver.BuildGlobalMatrices(stiffnessProvider);

        private void RebuildKs()
        {
            //TODO: This will rebuild all the stiffnesses of all subdomains, if even one subdomain has MaterialsModified = true.
            //      Optimize this, by passing a flag foreach subdomain to solver.BuildGlobalSubmatrices().

            bool mustRebuild = false;
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                if (subdomain.StiffnessModified)
                {
                    mustRebuild = true;
                    break;
                }
            }
            if (mustRebuild) ks = solver.BuildGlobalMatrices(stiffnessProvider);
            foreach (ISubdomain subdomain in model.Subdomains) subdomain.ResetMaterialsModifiedProperty();
        }

        private void BuildMs() => ms = solver.BuildGlobalMatrices(massProvider);

        private void BuildCs() => cs = solver.BuildGlobalMatrices(dampingProvider);

        private void BuildQs()
        {
            qs = new Dictionary<int, CsrMatrix>(model.SubdomainsDictionary.Count);
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                qs.Add(subdomain.ID, BuildQFromSubdomain(subdomain));
            }
        }

        //TODO: this should be done by an assembler class
        //TODO: make sure this is called whenever the ordering changes
        private CsrMatrix BuildQFromSubdomain(Subdomain subdomain) 
        {
            int numFreeDofs = subdomain.FreeDofOrdering.NumFreeDofs;
            var qSubdomain = DokRowMajor.CreateEmpty(numFreeDofs, numFreeDofs);
            DofTable allDofs = subdomain.FreeDofOrdering.FreeDofs;
            foreach (Element element in subdomain.Elements)
            {
                if (!(element.ElementType is IPorousFiniteElement)) continue;

                var e = (IPorousFiniteElement)element.ElementType;
                IMatrix q = e.CouplingMatrix(element);

                int iElementMatrixRow = 0;
                for (int i = 0; i < element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element).Count; i++)
                {
                    Node nodeRow = element.Nodes[i];
                    foreach (IDofType dofTypeRow in element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element)[i])
                    {
                        if (dofTypeRow != PorousMediaDof.Pressure) continue;

                        int dofRow = allDofs[nodeRow, dofTypeRow];
                        int iElementMatrixColumn = 0;

                        for (int j = 0; j < element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element).Count; j++)
                        {
                            Node nodeColumn = element.Nodes[j];
                            foreach (IDofType dofTypeColumn in element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element)[j])
                            {
                                if (dofTypeColumn == PorousMediaDof.Pressure) continue;

                                int dofColumn = allDofs[nodeColumn, dofTypeColumn];
                                qSubdomain.AddToEntry(dofColumn, dofRow, q[iElementMatrixRow, iElementMatrixColumn]);
                                iElementMatrixColumn++;
                            }
                        }
                        iElementMatrixRow++;
                    }
                }
            }

            return qSubdomain.BuildCsrMatrix(true);
        }

        private void ScaleSubdomainSolidVector(ISubdomain subdomain, IVector vector)
        {
            foreach ((INode node, IDofType dofType, int dofIdx) in subdomain.FreeDofOrdering.FreeDofs)
            {
                if (dofType != PorousMediaDof.Pressure) vector.Set(dofIdx, vector[dofIdx] * this.scalingCoefficient);
            }
        }

        private IMatrixView CalculateEffectiveMatrixInternal(ImplicitIntegrationCoefficients coefficients, ISubdomain subdomain)
        {
            int id = subdomain.ID;
            IMatrix matrix = this.Ks[id];
            matrix.LinearCombinationIntoThis(coefficients.Stiffness, Ms[id], coefficients.Mass);
            matrix.AxpyIntoThis(Cs[id], coefficients.Damping);
            return matrix;
        }

        #region IAnalyzerProvider Members
        public void ClearMatrices()
        {
            ms = null;
            cs = null;
            ks = null;
            qs = null;
        }

        public void Reset()
        {
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
                foreach (Element element in subdomain.Elements)
                    element.ElementType.ClearMaterialState();

            cs = null;
            ks = null;
            ms = null;
        }
        #endregion

        #region IImplicitIntegrationProvider Members

        public IMatrixView LinearCombinationOfMatricesIntoStiffness(ImplicitIntegrationCoefficients coefficients, 
            ISubdomain subdomain)
        {
            InitializeProvidersAndBuildQs(coefficients);
            scalingCoefficient = coefficients.Damping;
            return CalculateEffectiveMatrixInternal(coefficients, subdomain);
        }

        public void ProcessRhs(ImplicitIntegrationCoefficients coefficients, ISubdomain subdomain, IVector rhs)
        {
            scalingCoefficient = coefficients.Damping;
            ScaleSubdomainSolidVector(subdomain, rhs);
        }

        public IDictionary<int, IVector> GetRhsFromHistoryLoad(int timeStep)
        {
            foreach (Subdomain subdomain in model.Subdomains) subdomain.Forces.Clear(); //TODO: this is also done by model.AssignLoads()

            model.AssignLoads(solver.DistributeNodalLoads);
            model.AssignMassAccelerationHistoryLoads(timeStep);

            var rhsVectors = new Dictionary<int, IVector>();
            foreach (Subdomain subdomain in model.Subdomains) rhsVectors.Add(subdomain.ID, subdomain.Forces);
            return rhsVectors;
        }

        public IDictionary<int, IVector> GetAccelerationsOfTimeStep(int timeStep)
        {
            var d = new Dictionary<int, IVector>();
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                d.Add(linearSystem.Subdomain.ID, linearSystem.CreateZeroVector());
            }

            if (model.MassAccelerationHistoryLoads.Count > 0)
            {
                List<MassAccelerationLoad> m = new List<MassAccelerationLoad>(model.MassAccelerationHistoryLoads.Count);
                foreach (IMassAccelerationHistoryLoad l in model.MassAccelerationHistoryLoads)
                    m.Add(new MassAccelerationLoad() { Amount = l[timeStep], DOF = l.DOF });

                foreach (ISubdomain subdomain in model.Subdomains)
                {
                    int[] subdomainToGlobalDofs = model.GlobalDofOrdering.MapFreeDofsSubdomainToGlobal(subdomain);
                    foreach ((INode node, IDofType dofType, int subdomainDofIdx) in subdomain.FreeDofOrdering.FreeDofs)
                    {
                        int globalDofIdx = subdomainToGlobalDofs[subdomainDofIdx];
                        foreach (var l in m)
                        {
                            if (dofType == l.DOF) d[subdomain.ID].Set(globalDofIdx, l.Amount);
                        }
                    }
                }
            }

            //foreach (ElementMassAccelerationHistoryLoad load in model.ElementMassAccelerationHistoryLoads)
            //{
            //    MassAccelerationLoad hl = new MassAccelerationLoad() { Amount = load.HistoryLoad[timeStep] * 564000000, DOF = load.HistoryLoad.DOF };
            //    load.Element.Subdomain.AddLocalVectorToGlobal(load.Element,
            //        load.Element.ElementType.CalculateAccelerationForces(load.Element, (new MassAccelerationLoad[] { hl }).ToList()),
            //        load.Element.Subdomain.Forces);
            //}

            return d;
        }

        public IDictionary<int, IVector> GetVelocitiesOfTimeStep(int timeStep)
        {
            var d = new Dictionary<int, IVector>();

            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                d.Add(linearSystem.Subdomain.ID, linearSystem.CreateZeroVector());
            }

            return d;
        }

        public IVector MassMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
            => this.Ms[subdomain.ID].Multiply(vector);

        public IVector DampingMatrixVectorProduct(ISubdomain subdomain, IVectorView vector)
        {
            IVector result = this.Cs[subdomain.ID].Multiply(vector);
            result.SubtractIntoThis(qs[subdomain.ID].Multiply(vector));
            return result;
        }

        #endregion

        #region IStaticProvider Members

        public IMatrixView CalculateMatrix(ISubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        public (IMatrixView matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr) CalculateSubMatrices(ISubdomain subdomain)
        {
            throw new NotImplementedException();
        }

        //public void CalculateMatrices()
        //{
        //    throw new NotImplementedException();
        //}

        #endregion

        #region INonLinearProvider Members

        public double CalculateRhsNorm(IVectorView rhs)
        {
            //TODO: cache the relevant indices.
            //return (new Vector<double>(rhs)).Norm;

            double norm = 0;
            foreach ((INode node, IDofType dofType, int dofIdx) in model.GlobalDofOrdering.GlobalFreeDofs)
            {
                if (dofType != PorousMediaDof.Pressure) norm += rhs[dofIdx] * rhs[dofIdx];
            }
            return Math.Sqrt(norm);
        }

        public void ProcessInternalRhs(ISubdomain subdomain, IVectorView solution, IVector rhs)
        {
            rhs.AddIntoThis(qs[subdomain.ID].Multiply(solution));
            ScaleSubdomainSolidVector(subdomain, rhs);
        }

        #endregion
    }
}
