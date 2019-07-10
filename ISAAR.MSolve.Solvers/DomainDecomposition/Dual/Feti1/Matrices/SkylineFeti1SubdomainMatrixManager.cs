using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: Add state checking for all the managed matrices. A state machine (using State pattern) should help.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Matrices
{
    /// <summary>
    /// Dense format for Kbb, skyline for Kff, Kii and CSC for Kib.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SkylineFeti1SubdomainMatrixManager : IFeti1SubdomainMatrixManager
    {
        private readonly SkylineAssembler assembler = new SkylineAssembler();
        private readonly SingleSubdomainSystem<SkylineMatrix> linearSystem;
        private readonly IReorderingAlgorithm reordering;

        private SemidefiniteCholeskySkyline inverseKff;
        private LdlSkyline inverseKii;
        private DiagonalMatrix inverseKiiDiagonal;
        private Matrix Kbb;
        private CscMatrix Kib;

        public SkylineFeti1SubdomainMatrixManager(ISubdomain subdomain, IReorderingAlgorithm reordering)
        {
            this.linearSystem = new SingleSubdomainSystem<SkylineMatrix>(subdomain);
            this.reordering = reordering;
        }

        public ISingleSubdomainLinearSystem LinearSystem => linearSystem;

        public List<Vector> RigidBodyModes { get; private set; } = null;

        public IMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement> elements,
            IElementMatrixProvider matrixProvider)
            => assembler.BuildGlobalMatrix(dofOrdering, elements, matrixProvider);

        public (IMatrix Kff, IMatrixView Kfc, IMatrixView Kcf, IMatrixView Kcc) BuildGlobalSubmatrices(
            ISubdomainFreeDofOrdering freeDofOrdering, ISubdomainConstrainedDofOrdering constrainedDofOrdering,
            IEnumerable<IElement> elements, IElementMatrixProvider matrixProvider)
            => assembler.BuildGlobalSubmatrices(freeDofOrdering, constrainedDofOrdering, elements, matrixProvider);

        /// <summary>
        /// If the matrices stored in this object have already been calculated, they will be reused even if the original  
        /// free-free stiffness matrix has changed. To avoid that, this method must be called. 
        /// </summary>
        public void Clear()
        {
            inverseKff = null;
            RigidBodyModes = null;
            inverseKii = null;
            inverseKiiDiagonal = null;
            Kbb = null;
            Kib = null;
            //linearSystem.Matrix = null; // DO NOT DO THAT!!! The analyzer manages that.
        }

        /// <summary>
        /// Will do nothing if it was already called. To perform this for a different stiffness matrix, first call 
        /// <see cref="Clear"/>.
        /// </summary>
        public void ExtractAndInvertKiiDiagonal(int[] internalDofs)
        {
            if (inverseKiiDiagonal != null) return;
            try
            {
                var diagonal = new double[internalDofs.Length];
                for (int i = 0; i < diagonal.Length; ++i)
                {
                    int idx = internalDofs[i];
                    diagonal[i] = 1.0 / linearSystem.Matrix[idx, idx];
                    //diagonal[i] = linearSystem.Matrix[idx, idx];
                }
                inverseKiiDiagonal = DiagonalMatrix.CreateFromArray(diagonal, false);
                //inverseKiiDiagonal.Invert();
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The free-free stiffness matrix of this subdomain has been overwritten and cannot be used anymore."
                    + "Try calling this method before factorizing/inverting it.");
            }
        }

        /// <summary>
        /// Will do nothing if it was already called. To perform this for a different stiffness matrix, first call 
        /// <see cref="Clear"/>.
        /// </summary>
        public void ExtractAndInvertKii(int[] internalDofs)
        {
            if (inverseKii != null) return;
            try
            {
                SkylineMatrix Kii = linearSystem.Matrix.GetSubmatrixSymmetricSkyline(internalDofs);
                inverseKii = Kii.FactorLdl(true);
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The free-free stiffness matrix of this subdomain has been overwritten and cannot be used anymore."
                    + "Try calling this method before factorizing/inverting it.");
            }
        }

        /// <summary>
        /// Will do nothing if it was already called. To perform this for a different stiffness matrix, first call 
        /// <see cref="Clear"/>.
        /// </summary>
        public void ExtractKbb(int[] boundaryDofs)
        {
            if (Kbb != null) return;
            try
            {
                Kbb = linearSystem.Matrix.GetSubmatrixFull(boundaryDofs, boundaryDofs);
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The free-free stiffness matrix of this subdomain has been overwritten and cannot be used anymore."
                    + "Try calling this method before factorizing/inverting it.");
            }
        }

        /// <summary>
        /// Will do nothing if it was already called. To perform this for a different stiffness matrix, first call 
        /// <see cref="Clear"/>.
        /// </summary>
        public void ExtractKbiKib(int[] boundaryDofs, int[] internalDofs)
        {
            if (Kib != null) return;
            try
            {
                Kib = linearSystem.Matrix.GetSubmatrixCsc(internalDofs, boundaryDofs);
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The free-free stiffness matrix of this subdomain has been overwritten and cannot be used anymore."
                    + "Try calling this method before factorizing/inverting it.");
            }
        }

        public void HandleDofOrderingWillBeModified() => assembler.HandleDofOrderingWillBeModified();

        /// <summary>
        /// Will do nothing if it was already called. To perform this for a different stiffness matrix, first call 
        /// <see cref="Clear"/>.
        /// </summary>
        public void InvertKff(double factorizationTolerance, bool inPlace)
        {
            if (RigidBodyModes != null) return;
            inverseKff = linearSystem.Matrix.FactorSemidefiniteCholesky(inPlace, factorizationTolerance);
            RigidBodyModes = new List<Vector>();
            foreach (double[] rbm in inverseKff.NullSpaceBasis)
            {
                RigidBodyModes.Add(Vector.CreateFromArray(rbm, false));
            }
        }

        public Vector MultiplyInverseKffTimes(Vector vector)
        {
            if (inverseKff == null) throw new InvalidOperationException(
                "The free-free stiffness matrix 'Kff' of this subdomain must be factorized first.");
            return inverseKff.MultiplyGeneralizedInverseMatrixTimesVector(vector);
        }

        public Vector MultiplyInverseKiiDiagonalTimes(Vector vector)
        {
            if (inverseKiiDiagonal == null)
            {
                throw new InvalidOperationException("The inverse of the diagonal of the internal-internal stiffness submatrix"
                    + " 'inv(diag(Kii))' of this subdomain must be calculated first.");
            }
            return inverseKiiDiagonal * vector;
        }

        public Matrix MultiplyInverseKiiDiagonalTimes(Matrix matrix)
        {
            if (inverseKiiDiagonal == null)
            {
                throw new InvalidOperationException("The inverse of the diagonal of the internal-internal stiffness submatrix"
                    + " 'inv(diag(Kii))' of this subdomain must be calculated first.");
            }
            return inverseKiiDiagonal * matrix;
        }

        public Vector MultiplyInverseKiiTimes(Vector vector)
        {
            if (inverseKii == null)
            {
                throw new InvalidOperationException("The inverse of the internal-internal stiffness submatrix"
                    + " 'inv(Kii)' of this subdomain must be calculated first.");
            }
            return inverseKii.SolveLinearSystem(vector);
        }

        public Matrix MultiplyInverseKiiTimes(Matrix matrix)
        {
            if (inverseKii == null)
            {
                throw new InvalidOperationException("The inverse of the internal-internal stiffness submatrix"
                    + " 'inv(Kii)' of this subdomain must be calculated first.");
            }
            return inverseKii.SolveLinearSystems(matrix);
        }

        public Vector MultiplyKbbTimes(Vector vector)
        {
            if (Kbb == null)
            {
                throw new InvalidOperationException(
                    "The boundary-boundary stiffness submatrix 'Kbb' of this subdomain must be calculated first.");
            }
            return Kbb * vector;
        }

        public Matrix MultiplyKbbTimes(Matrix matrix)
        {
            if (Kbb == null)
            {
                throw new InvalidOperationException(
                    "The boundary-boundary stiffness submatrix 'Kbb' of this subdomain must be calculated first.");
            }
            return Kbb * matrix;
        }

        public Vector MultiplyKbiTimes(Vector vector)
        {
            if (Kib == null)
            {
                throw new InvalidOperationException(
                    "The boundary-internal stiffness submatrix 'Kbi' of this subdomain must be calculated first.");
            }
            return Kib.Multiply(vector, true);
        }

        public Matrix MultiplyKbiTimes(Matrix matrix)
        {
            if (Kib == null)
            {
                throw new InvalidOperationException(
                    "The boundary-internal stiffness submatrix 'Kbi' of this subdomain must be calculated first.");
            }
            return Kib.MultiplyRight(matrix, true);
        }

        public Vector MultiplyKibTimes(Vector vector)
        {
            if (Kib == null)
            {
                throw new InvalidOperationException(
                    "The internal-boundary stiffness submatrix 'Kib' of this subdomain must be calculated first.");
            }
            return Kib.Multiply(vector);
        }

        public Matrix MultiplyKibTimes(Matrix matrix)
        {
            if (Kib == null)
            {
                throw new InvalidOperationException(
                    "The internal-boundary stiffness submatrix 'Kib' of this subdomain must be calculated first.");
            }
            return Kib.MultiplyRight(matrix);
        }

        public void SetSolutionVector(Vector solution) => linearSystem.SolutionConcrete = solution;

        public void ReorderInternalDofs(Feti1DofSeparator dofSeparator, ISubdomain subdomain)
        {
            if (reordering == null) return; // Use the natural ordering and do not modify any stored dof data
            try
            {
                int[] internalDofs = dofSeparator.InternalDofIndices[subdomain.ID];
                var pattern = linearSystem.Matrix.GetSubmatrixSymmetricPattern(internalDofs);
                (int[] permutation, bool oldToNew) = reordering.FindPermutation(pattern);
                int[] newInternalDofs = ReorderingUtilities.ReorderKeysOfDofIndicesMap(internalDofs, permutation, oldToNew);

                // What if the dof separator gets added other state that needs to be updated?
                dofSeparator.InternalDofIndices[subdomain.ID] = newInternalDofs;
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The free-free stiffness matrix of this subdomain has been overwritten and cannot be used anymore."
                    + "Try calling this method before factorizing/inverting it.");
            }
        }

        public class Factory : IFeti1SubdomainMatrixManagerFactory
        {
            private readonly IReorderingAlgorithm reordering;

            //TODO: Use the reordering classes of project Solvers.
            //TODO: If the natural ordering is best, then there is no need to modify the stored dof data. 
            //      Find a better way to handle it, perhaps by checking if IReorderingAlgorithm produced a better ordering.
            public Factory(IReorderingAlgorithm reordering = null)
            {
                this.reordering = reordering;
            }

            public IFeti1SubdomainMatrixManager CreateMatricesManager(ISubdomain subdomain)
                => new SkylineFeti1SubdomainMatrixManager(subdomain, reordering);
        }
    }
}
