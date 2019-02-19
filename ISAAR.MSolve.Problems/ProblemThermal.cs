using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using System.Threading.Tasks;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using System;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.Solvers.Assemblers;

namespace ISAAR.MSolve.Problems
{
    public class ProblemThermal: IImplicitIntegrationProvider, IStaticProvider, INonLinearProvider
    {
        private Dictionary<int, IMatrix2D> capacityMatrices, conductivityMatrices;
        //private readonly IStructuralModel model; // We cannot generalize it yet.
        private readonly Model model;
        private IDictionary<int, ILinearSystem> subdomains;
        private ElementStructuralStiffnessProvider conductivityProvider = new ElementStructuralStiffnessProvider();
        private ElementStructuralMassProvider capacityProvider = new ElementStructuralMassProvider();

        public ProblemThermal(Model model, IDictionary<int, ILinearSystem> subdomains)
        {
            this.model = model;
            this.subdomains = subdomains;
        }

        private IDictionary<int, IMatrix2D> CapacityMatrices
        {
            get
            {
                if (capacityMatrices == null) BuildCapacityMatrices();
                return capacityMatrices;
            }
        }

        private IDictionary<int, IMatrix2D> ConductivityMatrices
        {
            get
            {
                if (conductivityMatrices == null) BuildConductivityMatrices();
                //else RebuildConductivityMatrices();
                return conductivityMatrices;
            }
        }

        private void BuildConductivityMatrices()
        {
            conductivityMatrices = new Dictionary<int, IMatrix2D>(model.ISubdomainsDictionary.Count);
            //ks.Add(1, new SkylineMatrix2D<double>(new double[,] { { 6, -2 }, { -2, 4 } }));
            ElementStructuralStiffnessProvider s = new ElementStructuralStiffnessProvider();
            //foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            //    ks.Add(subdomain.ID, GlobalMatrixAssemblerSkyline.CalculateGlobalMatrix(subdomain, s));

            //var kks = new Dictionary<int, IMatrix2D<double>>(model.SubdomainsDictionary.Count);
            int procs = VectorExtensions.AffinityCount;
            var k = model.ISubdomainsDictionary.Keys.Select(x => x).ToArray<int>();
            var internalKs = new Dictionary<int, IMatrix2D>[procs];
            Parallel.ForEach(k.PartitionLimits(procs), limit =>
            {
                if (limit.Item3 - limit.Item2 > 0)
                {
                    internalKs[limit.Item1] = new Dictionary<int, IMatrix2D>(limit.Item3 - limit.Item2);
                    for (int i = limit.Item2; i < limit.Item3; i++)
                        internalKs[limit.Item1].Add(k[i], GlobalMatrixAssemblerSkyline.CalculateFreeFreeGlobalMatrix(model.ISubdomainsDictionary[k[i]], s));
                }
                else
                    internalKs[limit.Item1] = new Dictionary<int, IMatrix2D>();
            });
            for (int i = 0; i < procs; i++)
                foreach (int key in internalKs[i].Keys)
                    conductivityMatrices.Add(key, internalKs[i][key]);
        }

        private void RebuildConductivityMatrices()
        {
            foreach (ISubdomain subdomain in model.ISubdomainsDictionary.Values)
            //Parallel.ForEach(model.SubdomainsDictionary.Values, subdomain =>
            {
                if (subdomain.MaterialsModified)
                {
                    conductivityMatrices[subdomain.ID] = GlobalMatrixAssemblerSkyline.CalculateFreeFreeGlobalMatrix(subdomain, conductivityProvider);
                    subdomain.ResetMaterialsModifiedProperty();
                }
            }
        }

        private void BuildCapacityMatrices()
        {
            capacityMatrices = new Dictionary<int, IMatrix2D>(model.ISubdomainsDictionary.Count);
            //ms.Add(1, new SkylineMatrix2D<double>(new double[,] { { 2, 0 }, { 0, 1 } }));
            ElementStructuralMassProvider s = new ElementStructuralMassProvider();
            foreach (ISubdomain subdomain in model.ISubdomainsDictionary.Values)
                capacityMatrices.Add(subdomain.ID, GlobalMatrixAssemblerSkyline.CalculateFreeFreeGlobalMatrix(subdomain, s));
        }

        #region IAnalyzerProvider Members
        public void Reset()
        {
            foreach (Subdomain subdomain in model.ISubdomainsDictionary.Values)
                foreach (var element in subdomain.ElementsDictionary.Values)
                    element.ElementType.ClearMaterialState();

            conductivityMatrices = null;
            capacityMatrices = null;
        }
        #endregion 

        #region IImplicitIntegrationProvider Members

        public void CalculateEffectiveMatrix(ILinearSystem subdomain, ImplicitIntegrationCoefficients coefficients)
        {
            // Either way the GlobalMatrixAssemblerSkyline is hardcoded at this point, so the matrix, can only be Skyline.
            SkylineMatrix2D conductivity = (SkylineMatrix2D)(this.ConductivityMatrices[subdomain.ID]);

            // We will keep both the conductivity and the effective matrix, therefore only the effective matrix must be factorized. 
            // CalculateMatrix() set the conductivity as the linear system matrix and then the solver will factorize it.
            //if (((SkylineMatrix2D)subdomain.Matrix).IsFactorized) BuildConductivityMatrices();
            if (conductivity.IsFactorized) throw new InvalidOperationException("Conductivity matrix has been factorized. This should not have happened in a dynamic setting.");

            SkylineMatrix2D effective = conductivity.Copy(); // We need the conductivity matrix in the analyzer if the integration scheme is not purely implicit.
            effective.LinearCombination(
                new double[]
                {
                    coefficients.Stiffness, coefficients.Mass
                },
                new IMatrix2D[]
                {
                    effective, this.CapacityMatrices[subdomain.ID]
                });

            subdomain.Matrix = effective;
        }

        public void ProcessRHS(ILinearSystem subdomain, ImplicitIntegrationCoefficients coefficients)
        {
            // Method intentionally left empty.
        }

        public IDictionary<int, double[]> GetAccelerationsOfTimeStep(int timeStep)
        {
            throw new InvalidOperationException("This is does not make sense in explicit methods for first order equations");
        }

        //TODO: not used in explicit methods for first order equations
        public IDictionary<int, double[]> GetVelocitiesOfTimeStep(int timeStep)
        {
            throw new InvalidOperationException("This is not needed in explicit methods for first order equations");
        }

        //TODO: It is very common for the loading to be constant w.r.t. time throughout the analysis. This results in the same
        //      rhs at each iteration and we could avoid many operations. Find a way to allow this optimization.
        public void GetRHSFromHistoryLoad(int timeStep)
        {
            foreach (Subdomain subdomain in model.ISubdomainsDictionary.Values)
            {
                Array.Clear(subdomain.Forces, 0, subdomain.Forces.Length);
            }

            model.AssignNodalLoads(); // Time-independent nodal loads
            //model.AssignTimeDependentNodalLoads(timeStep); // Time-dependent nodal loads

            foreach (var l in subdomains)
                l.Value.RHS.CopyFrom(0, l.Value.RHS.Length, new Vector(model.ISubdomainsDictionary[l.Key].Forces), 0);

        }

        public void MassMatrixVectorProduct(ILinearSystem subdomain, IVector vIn, double[] vOut)
        {
            this.capacityMatrices[subdomain.ID].Multiply(vIn, vOut);
        }

        public void DampingMatrixVectorProduct(ILinearSystem subdomain, IVector vIn, double[] vOut)
        {
            this.conductivityMatrices[subdomain.ID].Multiply(vIn, vOut);
        }

        #endregion

        #region IStaticProvider Members

        public void CalculateMatrix(ILinearSystem subdomain)
        {
            if (conductivityMatrices == null) BuildConductivityMatrices();
            subdomain.Matrix = this.conductivityMatrices[subdomain.ID];
        }

        #endregion

        #region INonLinearProvider Members

        public double RHSNorm(double[] rhs)
        {
            return (new Vector(rhs)).Norm;
        }

        public void ProcessInternalRHS(ILinearSystem subdomain, double[] rhs, double[] solution)
        {
        }

        #endregion
    }
}
