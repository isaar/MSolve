using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using System.Linq;

//TODO: extend it to 3D.
//TODO: What about passive elements? They should not be included in the analysis.
//TODO: Add support for various boundary conditions and load cases
//TODO: When the general classes achieve the efficiency of the 2D uniform hardcoded classes (such as this), then the latter
//      must be removed.
namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Analysis
{
    /// <summary>
    /// This class implements the FEM analysis part of the 99 line topology code, which is 2D, linear elastic with a uniform mesh.
    /// In this case many optimizations are possible (e.g. only calculating the element stiffness matrix once.). Therefore this 
    /// class is useful to compare other FEM analyses with, in terms of computational cost. 
    /// For the full paper see "A 99 line topology optimization code written in Matlab, O. Sigmund, 1991"
    /// </summary>
    public class LinearFemAnalysis2DUniformHardcoded: ILinearFemAnalysis
    {
        //TODO: Use MSolve analysis and Model instead of this
        public enum BoundaryConditions
        {
            MbbBeam, ShortCantilever, Cantilever2LoadCases
        }
        private readonly BoundaryConditions bc;

        private const int numElementDofs = 8;
        private static readonly int[] elementDofsLocal = { 0, 1, 2, 3, 4, 5, 6, 7 };

        private readonly int numElementsX, numElementsY;
        private readonly ElasticMaterial2D material;

        private Vector[] globalDisplacements;
        private IMatrixView unitStiffness; //TODO: Should this be readonly without an Initialize() method?
        private IVectorView youngModuli;

        public LinearFemAnalysis2DUniformHardcoded(int numElementsX, int numElementsY, ElasticMaterial2D material,
            BoundaryConditions bc)
        {
            this.numElementsX = numElementsX;
            this.numElementsY = numElementsY;
            this.NumElements = numElementsX * numElementsY;
            this.material = material;
            this.bc = bc;
        }

        public int NumElements { get; }

        public int NumLoadCases => (bc == BoundaryConditions.Cantilever2LoadCases) ? 2 : 1;

        public void UpdateModelAndAnalyze(IVectorView youngModuli)
        {
            this.youngModuli = youngModuli;

            // Global stiffness matrix assembly
            int numAllDofs = 2 * (numElementsX + 1) * (numElementsY + 1);
            var K = DokSymmetric.CreateEmpty(numAllDofs);
            for (int e = 0; e < NumElements; ++e)
            {
                int[] elementDofsGlobal = GetElementDofs(e);
                K.AddSubmatrixSymmetric(unitStiffness.Scale(youngModuli[e]), elementDofsLocal, elementDofsGlobal);
            }

            // Apply boundary conditions
            (int[] freeDofs, Vector[] Fs) = ApplyBoundaryConditions();
            int numLoadCases = Fs.Length;
            globalDisplacements = new Vector[numLoadCases];

            // Solve linear system
            DokSymmetric Kf = K.GetSubmatrix(freeDofs);
            var factor = CholeskyCSparseNet.Factorize(Kf.BuildSymmetricCscMatrix(true));
            for (int c = 0; c < numLoadCases; ++c)
            {
                Vector Ff = Fs[c].GetSubvector(freeDofs);
                Vector Uf = factor.SolveLinearSystem(Ff);
                var U = Vector.CreateZero(numAllDofs);
                for (int i = 0; i < freeDofs.Length; ++i) U[freeDofs[i]] = Uf[i];
                globalDisplacements[c] = U;
                //U.CopyNonContiguouslyFrom(freeDofs, Uf, Enumerable.Range(0, freeDofs.Length).ToArray()); // alternative way, but probably slower.
            }
        }

        //TODO: doesn't this depend on the geometry? Or is it normalized somehow?
        public double CalculateTotalVolume(IVectorView densities) => densities.Sum();

        public IMatrixView GetElementStiffnessForCurrentMaterial(int elementIdx) 
            => unitStiffness.Scale(youngModuli[elementIdx]);

        public IMatrixView GetElementStiffnessForUnitMaterial(int elementIdx) => unitStiffness;

        public Vector GetElementDisplacements(int elementIdx, int loadCaseIdx)
        {
            if (globalDisplacements == null) throw new InvalidOperationException(
                "The global displacements must be calculated first");
            return globalDisplacements[loadCaseIdx].GetSubvector(GetElementDofs(elementIdx));
        }

        public void Initialize()
        {
            unitStiffness = CalcQuad4Stiffness();
        }

        private (int[] freeDofs, Vector[] forcesPerLoadCase) ApplyBoundaryConditions()
        {
            int numAllDofs = 2 * (numElementsX + 1) * (numElementsY + 1);
            if (bc == BoundaryConditions.MbbBeam)
            {
                // Define loads and supports for half MBB beam.
                var Fs = new Vector[] { Vector.CreateZero(numAllDofs) };
                Fs[0][1] = -1.0;
                var fixedDofs = new HashSet<int>(); //TODO: Use LINQ to simplify this
                for (int i = 0; i < 2 * (numElementsY + 1); i += 2) fixedDofs.Add(i);
                fixedDofs.Add(numAllDofs - 1);
                int[] freeDofs = Enumerable.Range(0, numAllDofs).Except(fixedDofs).ToArray();
                return (freeDofs, Fs);
            }
            else if (bc == BoundaryConditions.ShortCantilever)
            {
                // Define loads and supports for cantilever beam
                var Fs = new Vector[] { Vector.CreateZero(numAllDofs) };
                Fs[0][numAllDofs - 1] = -1.0;
                IEnumerable<int> fixedDofs = Enumerable.Range(0, 2 * (numElementsY + 1));
                int[] freeDofs = Enumerable.Range(0, numAllDofs).Except(fixedDofs).ToArray();
                return (freeDofs, Fs);
            }
            else if (bc == BoundaryConditions.Cantilever2LoadCases)
            {
                // Define supports for cantilever beam
                IEnumerable<int> fixedDofs = Enumerable.Range(0, 2 * (numElementsY + 1));
                int[] freeDofs = Enumerable.Range(0, numAllDofs).Except(fixedDofs).ToArray();

                // Load case 1: unit load towards -y at bottom right corner
                var Fs = new Vector[] { Vector.CreateZero(numAllDofs), Vector.CreateZero(numAllDofs) };
                Fs[0][2 * (numElementsX + 1) * (numElementsY + 1) - 1] = -1.0;

                // Load case 2: unit load towards +y at top right corner
                Fs[1][2 * (numElementsX) * (numElementsY + 1) + 1] = 1.0;
                return (freeDofs, Fs);
            }
            else throw new Exception("This code should not have been reached.");
        }

        private Matrix CalcQuad4Stiffness()
        {
            double E = material.YoungModulus;
            double nu = material.PoissonRatio;
            double[] k = { 0.5 - nu / 6.0, 0.125 + nu / 8.0, -0.25 - nu / 12.0, -0.125 + 3 * nu / 8.0,
                -0.25 + nu / 12.0, -0.125 - nu / 8.0, nu / 6.0, 0.125 - 3 * nu / 8.0 }; // unique stiffness matrix entries
            var Ke = Matrix.CreateFromArray(new double[,]
            {
                 { k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7] },
                 { k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2] },
                 { k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1] },
                 { k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4] },
                 { k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3] },
                 { k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6] },
                 { k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5] },
                 { k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0] }
            });
            Ke.ScaleIntoThis(E / (1 - nu * nu));
            return Ke;
        }

        /// <summary>
        /// The global indices of the element's dofs, in 0-based numbering.
        /// </summary>
        /// <param name="elx">0-based numbering</param>
        /// <param name="ely">0-based numbering</param>
        private int[] GetElementDofs(int elementIdx)
        {
            // 1D indexing (0-based) to 2D indexing (1-based, y-major)
            int ely = elementIdx / numElementsX + 1; 
            int elx = elementIdx % numElementsX + 1;

            int n1 = (numElementsY + 1) * (elx - 1) + ely;
            int n2 = (numElementsY + 1) * elx + ely;
            return new int[]
            {
                2 * n1 - 2,     2 * n1 - 1,     2 * n2 - 2,     2 * n2 - 1,
                2 * n2,         2 * n2 + 1,     2 * n1,         2 * n1 + 1
            };
        }
    }
}
