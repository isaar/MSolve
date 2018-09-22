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
    public class ProblemStructural : IImplicitIntegrationProvider, IStaticProvider, INonLinearProvider
    {
        private Dictionary<int, IMatrix2D> ms, cs, ks;
        private readonly IStructuralModel model;
        private IDictionary<int, ILinearSystem> subdomains;
        private ElementStructuralStiffnessProvider stiffnessProvider = new ElementStructuralStiffnessProvider();
        private ElementStructuralMassProvider massProvider = new ElementStructuralMassProvider();

        public ProblemStructural(IStructuralModel model, IDictionary<int, ILinearSystem> subdomains)
        {
            this.model = model;
            this.subdomains = subdomains;
        }

        public double AboserberE { get; set; }
        public double Aboseberv { get; set; }

        private IDictionary<int, IMatrix2D> Ms
        {
            get
            {
                if (ms == null) BuildMs();
                return ms;
            }
        }

        private IDictionary<int, IMatrix2D> Cs
        {
            get
            {
                if (cs == null) BuildCs();
                return cs;
            }
        }

        private IDictionary<int, IMatrix2D> Ks
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

        private void BuildKs()
        {
            ks = new Dictionary<int, IMatrix2D>(model.ISubdomainsDictionary.Count);
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
                        internalKs[limit.Item1].Add(k[i], GlobalMatrixAssemblerSkyline.CalculateGlobalMatrix(model.ISubdomainsDictionary[k[i]], s));
                }
                else
                    internalKs[limit.Item1] = new Dictionary<int, IMatrix2D>();
            });
            for (int i = 0; i < procs; i++)
                foreach (int key in internalKs[i].Keys)
                    ks.Add(key, internalKs[i][key]);
        }

        private void RebuildKs()
        {
            foreach (ISubdomain subdomain in model.ISubdomainsDictionary.Values)
            //Parallel.ForEach(model.SubdomainsDictionary.Values, subdomain =>
            {
                if (subdomain.MaterialsModified)
                {
                    ks[subdomain.ID] = GlobalMatrixAssemblerSkyline.CalculateGlobalMatrix(subdomain, stiffnessProvider);
                    subdomain.ResetMaterialsModifiedProperty();
                }
            }
        }

        private void BuildMs()
        {
            ms = new Dictionary<int, IMatrix2D>(model.ISubdomainsDictionary.Count);
            //ms.Add(1, new SkylineMatrix2D<double>(new double[,] { { 2, 0 }, { 0, 1 } }));
            ElementStructuralMassProvider s = new ElementStructuralMassProvider();
            foreach (ISubdomain subdomain in model.ISubdomainsDictionary.Values)
                ms.Add(subdomain.ID, GlobalMatrixAssemblerSkyline.CalculateGlobalMatrix(subdomain, s));
        }

        private void BuildCs()
        {
            cs = new Dictionary<int, IMatrix2D>(model.ISubdomainsDictionary.Count);
            //foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            //    cs.Add(subdomain.ID, SkylineMatrix2D<double>.Empty(subdomain.TotalDOFs));
            ElementStructuralDampingProvider s = new ElementStructuralDampingProvider();
            foreach (ISubdomain subdomain in model.ISubdomainsDictionary.Values)
                cs.Add(subdomain.ID, GlobalMatrixAssemblerSkyline.CalculateGlobalMatrix(subdomain, s));
        }

        #region IAnalyzerProvider Members
        public void Reset()
        {
            foreach (Subdomain subdomain in model.ISubdomainsDictionary.Values)
                foreach (var element in subdomain.ElementsDictionary.Values)
                    element.ElementType.ClearMaterialState();

            cs = null;
            ks = null;
            ms = null;
        }
        #endregion 

        #region IImplicitIntegrationProvider Members

        public void CalculateEffectiveMatrix(ILinearSystem subdomain, ImplicitIntegrationCoefficients coefficients)
        {
            subdomain.Matrix = this.Ks[subdomain.ID];

            if (((SkylineMatrix2D)subdomain.Matrix).IsFactorized)
                BuildKs();

            var m = subdomain.Matrix as ILinearlyCombinable;
            m.LinearCombination(
                new double[]
                {
                    coefficients.Stiffness, coefficients.Mass, coefficients.Damping
                },
                new IMatrix2D[]
                {
                    this.Ks[subdomain.ID], this.Ms[subdomain.ID], this.Cs[subdomain.ID]
                });
        }

        public void ProcessRHS(ILinearSystem subdomain, ImplicitIntegrationCoefficients coefficients)
        {
            // Method intentionally left empty.
        }

        public IDictionary<int, double[]> GetAccelerationsOfTimeStep(int timeStep)
        {
            var d = new Dictionary<int, double[]>();
            foreach (Subdomain subdomain in model.ISubdomainsDictionary.Values)
                d.Add(subdomain.ID, new double[subdomain.TotalDOFs]);

            if (model.MassAccelerationHistoryLoads.Count > 0)
            {
                List<MassAccelerationLoad> m = new List<MassAccelerationLoad>(model.MassAccelerationHistoryLoads.Count);
                foreach (IMassAccelerationHistoryLoad l in model.MassAccelerationHistoryLoads)
                    m.Add(new MassAccelerationLoad() { Amount = l[timeStep], DOF = l.DOF });

                foreach (ISubdomain subdomain in model.ISubdomainsDictionary.Values)
                {
                    foreach (var nodeInfo in subdomain.GlobalNodalDOFsDictionary)
                    {
                        foreach (var dofPair in nodeInfo.Value)
                        {
                            foreach (var l in m)
                            {
                                if (dofPair.Key == l.DOF&& dofPair.Value!=-1)
                                {
                                    d[subdomain.ID][dofPair.Value] = l.Amount;
                                }
                            }
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

        public IDictionary<int, double[]> GetVelocitiesOfTimeStep(int timeStep)
        {
            var d = new Dictionary<int, double[]>();
            foreach (Subdomain subdomain in model.ISubdomainsDictionary.Values)
                d.Add(subdomain.ID, new double[subdomain.TotalDOFs]);

            return d;
        }

        public void GetRHSFromHistoryLoad(int timeStep)
        {


            foreach (Subdomain subdomain in model.ISubdomainsDictionary.Values)
                for (int i = 0; i < subdomain.Forces.Length; i++)
                    subdomain.Forces[i] = 0;


            model.AssignLoads();
            model.AssignMassAccelerationHistoryLoads(timeStep);

            foreach (var l in subdomains)
                l.Value.RHS.CopyFrom(0, l.Value.RHS.Length, new Vector(model.ISubdomainsDictionary[l.Key].Forces), 0);

            ////AMBROSIOS
            //if (model.MassAccelerationHistoryLoads.Count > 0)
            //{
            //    int[] sDOFs = new int[] { 0, 1, 2, 3, 4, 5, 36, 37, 38, 39, 40, 41 };
            //    int[] fDOFs = new int[] { model.NodalDOFsDictionary[3][DOFType.X], model.NodalDOFsDictionary[3][DOFType.Y], model.NodalDOFsDictionary[3][DOFType.Z], 
            //        model.NodalDOFsDictionary[4][DOFType.X], model.NodalDOFsDictionary[4][DOFType.Y], model.NodalDOFsDictionary[4][DOFType.Z], 
            //        model.NodalDOFsDictionary[15][DOFType.X], model.NodalDOFsDictionary[15][DOFType.Y], model.NodalDOFsDictionary[16][DOFType.Z], 
            //        model.NodalDOFsDictionary[16][DOFType.X], model.NodalDOFsDictionary[15][DOFType.Y], model.NodalDOFsDictionary[16][DOFType.Z]};
            //    int[] fDOFs = new int[] { 6, 7, 8, 9, 10, 11, 42, 43, 44, 45, 46, 47 };
            //    var msfdata = new double[sDOFs.Length, fDOFs.Length];
            //    var lines = File.ReadAllLines(@"C:\Development\vs 2010\AnalSharp\msfAmbrosios.csv");
            //    for (int i = 0; i < sDOFs.Length; i++)
            //        for (int j = 0; j < fDOFs.Length; j++)
            //            msfdata[i, j] = Double.Parse(lines[i * fDOFs.Length + j]); 
            //            msfdata[i, j] = ms[1][sDOFs[i], fDOFs[j]];
            //    var msf = new Matrix2D<double>(msfdata);
            //    var acc = new double[sDOFs.Length];

            //    List<MassAccelerationLoad> m = new List<MassAccelerationLoad>(model.MassAccelerationHistoryLoads.Count);
            //    foreach (IMassAccelerationHistoryLoad l in model.MassAccelerationHistoryLoads)
            //        for (int i = 0; i < sDOFs.Length / 3; i++)
            //            acc[3*i] = l[timeStep];

            //    var localForces = new double[sDOFs.Length];
            //    msf.Multiply(new Vector<double>(acc), localForces);

            //    for (int j = 0; j < fDOFs.Length; j++)
            //        model.SubdomainsDictionary[1].Forces[fDOFs[j]] -= localForces[j];
            //}
        }

        public void MassMatrixVectorProduct(ILinearSystem subdomain, IVector vIn, double[] vOut)
        {
            this.Ms[subdomain.ID].Multiply(vIn, vOut);
        }

        public void DampingMatrixVectorProduct(ILinearSystem subdomain, IVector vIn, double[] vOut)
        {
            this.Cs[subdomain.ID].Multiply(vIn, vOut);
        }

        #endregion

        #region IStaticProvider Members

        public void CalculateMatrix(ILinearSystem subdomain)
        {
            if (ks == null) BuildKs();
            subdomain.Matrix = this.ks[subdomain.ID];
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
