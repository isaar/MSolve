using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.SchurComplements;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

//TODO: Kff should probably be a DOK. It will only be used to extract Krr, Krc, Kcc. 
//      What about dynamic problems, where Kff needs to do linear combinations and matrix-vector multiplications
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices
{
    /// <summary>
    /// Dense format for Kbb, Kcc, KccStar, skyline for Kff, Krr, Kii and CSC for Kib, Krc.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SkylineFetiDPSubdomainMatrixManager : IFetiDPSubdomainMatrixManager
    {
        private readonly SkylineAssembler assembler = new SkylineAssembler();
        private readonly SingleSubdomainSystem<SkylineMatrix> linearSystem;
        private readonly IReorderingAlgorithm reordering;

        private DiagonalMatrix inverseKiiDiagonal;
        private LdlSkyline inverseKii;
        private LdlSkyline inverseKrr;
        private Matrix Kbb;
        private CscMatrix Kib;
        private Matrix Kcc;
        private Matrix KccStar;
        private CscMatrix Krc;
        private SkylineMatrix Krr;

        public SkylineFetiDPSubdomainMatrixManager(ISubdomain subdomain, IReorderingAlgorithm reordering)
        {
            this.linearSystem = new SingleSubdomainSystem<SkylineMatrix>(subdomain);
            this.reordering = reordering;
        }

        public ISingleSubdomainLinearSystem LinearSystem => linearSystem;

        public Matrix SchurComplementOfRemainderDofs => KccStar;

        public IMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement> elements,
            IElementMatrixProvider matrixProvider)
            => assembler.BuildGlobalMatrix(dofOrdering, elements, matrixProvider);


        public (IMatrix Kff, IMatrixView Kfc, IMatrixView Kcf, IMatrixView Kcc) BuildGlobalSubmatrices(
            ISubdomainFreeDofOrdering freeDofOrdering, ISubdomainConstrainedDofOrdering constrainedDofOrdering,
            IEnumerable<IElement> elements, IElementMatrixProvider matrixProvider)
            => assembler.BuildGlobalSubmatrices(freeDofOrdering, constrainedDofOrdering, elements, matrixProvider);

        /// <summary>
        /// Will do nothing if it was already called. To perform this for a different stiffness matrix, first call 
        /// <see cref="Clear"/>.
        /// </summary>
        public void Clear()
        {
            inverseKii = null;
            inverseKiiDiagonal = null;
            inverseKrr = null;
            Kbb = null;
            Kib = null;
            Kcc = null;
            Krc = null;
            Krr = null;
            KccStar = null;
            //linearSystem.Matrix = null; // DO NOT DO THAT!!! The analyzer manages that.
        }

        /// <summary>
        /// Will do nothing if it was already called. To perform this for a different stiffness matrix, first call 
        /// <see cref="Clear"/>.
        /// </summary>
        public void CalcSchurComplementOfRemainderDofs()
        {
            // KccStar[s] = Kcc[s] - Krc[s]^T * inv(Krr[s]) * Krc[s]
            if (KccStar != null) return;
            KccStar = SchurComplementCsc.CalcSchurComplementFull(Kcc, Krc, inverseKrr);
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
                SkylineMatrix Kii = Krr.GetSubmatrixSymmetricSkyline(internalDofs);
                inverseKii = Kii.FactorLdl(true);
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The remainder-remainder stiffness submatrix of this subdomain has been overwritten and cannot be used"
                    + " anymore. Try calling this method before factorizing/inverting it.");
            }
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
                    diagonal[i] = 1.0 / Krr[idx, idx];
                    //diagonal[i] = Krr[idx, idx];
                }
                inverseKiiDiagonal = DiagonalMatrix.CreateFromArray(diagonal, false);
                //inverseKiiDiagonal.Invert();
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The remainder-remainder stiffness submatrix of this subdomain has been overwritten and cannot be used"
                    + " anymore. Try calling this method before factorizing/inverting it.");
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
                Kbb = Krr.GetSubmatrixSymmetricFull(boundaryDofs);
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The remainder-remainder stiffness submatrix of this subdomain has been overwritten and cannot be used"
                    + " anymore. Try calling this method before factorizing/inverting it.");
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
                Kib = Krr.GetSubmatrixCsc(internalDofs, boundaryDofs);
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The remainder-remainder stiffness submatrix of this subdomain has been overwritten and cannot be used"
                    + " anymore. Try calling this method before factorizing/inverting it.");
            }
        }

        /// <summary>
        /// Will do nothing if it was already called. To perform this for a different stiffness matrix, first call 
        /// <see cref="Clear"/>.
        /// </summary>
        public void ExtractKcc(int[] cornerDofs)
        {
            if (Kcc != null) return;
            Kcc = linearSystem.Matrix.GetSubmatrixFull(cornerDofs, cornerDofs);
        }

        /// <summary>
        /// Will do nothing if it was already called. To perform this for a different stiffness matrix, first call 
        /// <see cref="Clear"/>.
        /// </summary>
        public void ExtractKcrKrc(int[] cornerDofs, int[] remainderDofs)
        {
            if (Krc != null) return;
            Krc = linearSystem.Matrix.GetSubmatrixCsc(remainderDofs, cornerDofs);
        }

        /// <summary>
        /// Will do nothing if it was already called. To perform this for a different stiffness matrix, first call 
        /// <see cref="Clear"/>.
        /// </summary>
        public void ExtractKrr(int[] remainderDofs)
        {
            if (Krr != null)
            {
                if (inverseKrr != null)
                {
                    throw new InvalidOperationException("The remainder-remainder stiffness submatrix of this subdomain has"
                        + " already been calculated and then overwritten and cannot be used anymore. Restructure your code so"
                        + " that all operations that need Krr are finished before inverting it.");
                }
                return;
            }
            Krr = linearSystem.Matrix.GetSubmatrixSymmetricSkyline(remainderDofs);
        }

        public void HandleDofOrderingWillBeModified() => assembler.HandleDofOrderingWillBeModified();

        /// <summary>
        /// Will do nothing if it was already called. To perform this for a different stiffness matrix, first call 
        /// <see cref="Clear"/>.
        /// </summary>
        public void InvertKrr(bool inPlace)
        {
            if (inverseKrr != null) return;
            inverseKrr = Krr.FactorLdl(inPlace);
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
                throw new InvalidOperationException("The inverse of the internal_remainder - internal_remainder stiffness"
                    + " submatrix 'inv(Kii)' of this subdomain must be calculated first.");
            }
            return inverseKii.SolveLinearSystem(vector);
        }

        public Matrix MultiplyInverseKiiTimes(Matrix matrix)
        {
            if (inverseKii == null)
            {
                throw new InvalidOperationException("The inverse of the internal_remainder - internal_remainder stiffness"
                    + " submatrix 'inv(Kii)' of this subdomain must be calculated first.");
            }
            return inverseKii.SolveLinearSystems(matrix);
        }

        public Vector MultiplyInverseKrrTimes(Vector vector)
        {
            if (inverseKrr == null) throw new InvalidOperationException("The inverse of the remainder-remainder stiffness "
                + " submatrix 'Krr' of this subdomain must be calculated first.");
            return inverseKrr.SolveLinearSystem(vector);
        }

        public Vector MultiplyKbbTimes(Vector vector)
        {
            if (Kbb == null)
            {
                throw new InvalidOperationException("The boundary_remainder - boundary_remainder stiffness submatrix"
                    + " 'Kbb' of this subdomain must be calculated first.");
            }
            return Kbb * vector;
        }

        public Matrix MultiplyKbbTimes(Matrix matrix)
        {
            if (Kbb == null)
            {
                throw new InvalidOperationException("The boundary_remainder - boundary_remainder stiffness submatrix"
                    + " 'Kbb' of this subdomain must be calculated first.");
            }
            return Kbb * matrix;
        }

        public Vector MultiplyKbiTimes(Vector vector)
        {
            if (Kib == null)
            {
                throw new InvalidOperationException("The boundary_remainder - internal_remainder stiffness submatrix"
                    + " 'Kbi' of this subdomain must be calculated first.");
            }
            return Kib.Multiply(vector, true);
        }

        public Matrix MultiplyKbiTimes(Matrix matrix)
        {
            if (Kib == null)
            {
                throw new InvalidOperationException("The boundary_remainder - internal_remainder stiffness submatrix"
                    + " 'Kbi' of this subdomain must be calculated first.");
            }
            return Kib.MultiplyRight(matrix, true);
        }

        public Vector MultiplyKccTimes(Vector vector)
        {
            if (Kcc == null)
            {
                throw new InvalidOperationException("The corner-corner stiffness submatrix"
                    + " 'Kcc' of this subdomain must be calculated first.");
            }
            return Kcc * vector;
        }

        public Vector MultiplyKcrTimes(Vector vector)
        {
            if (Krc == null)
            {
                throw new InvalidOperationException("The remainder-corner stiffness submatrix"
                    + " 'Kcr' of this subdomain must be calculated first.");
            }
            return Krc.Multiply(vector, true);
        }

        public Vector MultiplyKibTimes(Vector vector)
        {
            if (Kib == null)
            {
                throw new InvalidOperationException("The internal_remainder - boundary_remainder stiffness submatrix"
                    + " 'Kib' of this subdomain must be calculated first.");
            }
            return Kib.Multiply(vector);
        }

        public Matrix MultiplyKibTimes(Matrix matrix)
        {
            if (Kib == null)
            {
                throw new InvalidOperationException("The internal_remainder - boundary_remainder stiffness submatrix"
                    + " 'Kib' of this subdomain must be calculated first.");
            }
            return Kib.MultiplyRight(matrix);
        }

        public Vector MultiplyKrcTimes(Vector vector)
        {
            if (Krc == null)
            {
                throw new InvalidOperationException("The remainder-corner stiffness submatrix"
                    + " 'Krc' of this subdomain must be calculated first.");
            }
            return Krc.Multiply(vector);
        }

        public void ReorderInternalDofs(FetiDPDofSeparator dofSeparator, ISubdomain subdomain)
        {
            if (reordering == null) return; // Use the natural ordering and do not modify any stored dof data
            try
            {
                int[] internalDofs = dofSeparator.InternalDofIndices[subdomain.ID];
                var pattern = Krr.GetSubmatrixSymmetricPattern(internalDofs);
                (int[] permutation, bool oldToNew) = reordering.FindPermutation(pattern);
                int[] newInternalDofs = ReorderingUtilities.ReorderKeysOfDofIndicesMap(internalDofs, permutation, oldToNew);

                // What if the dof separator gets added other state that needs to be updated?
                dofSeparator.InternalDofIndices[subdomain.ID] = newInternalDofs;
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The remainder-remainder stiffness submatrix of this subdomain has been already been calculated and"
                    + " then overwritten and cannot be used anymore. Try calling this method before"
                    + " factorizing/inverting it.");
            }
        }

        public void ReorderRemainderDofs(FetiDPDofSeparator dofSeparator, ISubdomain subdomain)
        {
            if (reordering == null) return; // Use the natural ordering and do not modify any stored dof data
            try
            {
                int s = subdomain.ID;
                int[] remainderDofs = dofSeparator.RemainderDofIndices[s];
                var pattern = linearSystem.Matrix.GetSubmatrixSymmetricPattern(remainderDofs);
                (int[] permutation, bool oldToNew) = reordering.FindPermutation(pattern);
                int[] newRemainderDofs = ReorderingUtilities.ReorderKeysOfDofIndicesMap(remainderDofs, permutation, oldToNew);

                // What if the dof separator gets added other state that needs to be updated?
                dofSeparator.RemainderDofIndices[s] = newRemainderDofs;
                dofSeparator.RemainderDofOrderings[s].Reorder(permutation, oldToNew);
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The free-free stiffness matrix of this subdomain has been overwritten and cannot be used anymore."
                    + "Try calling this method before factorizing/inverting it.");
            }
        }

        public class Factory : IFetiDPSubdomainMatrixManagerFactory
        {
            private readonly IReorderingAlgorithm reordering;

            //TODO: Use the reordering classes of project Solvers.
            //TODO: If the natural ordering is best, then there is no need to modify the stored dof data. 
            //      Find a better way to handle it, perhaps by checking if IReorderingAlgorithm produced a better ordering.
            public Factory(IReorderingAlgorithm reordering = null)
            {
                this.reordering = reordering;
            }

            public IFetiDPCoarseProblemSolver CreateCoarseProblemSolver(IReadOnlyList<ISubdomain> subdomains)
                => new SkylineFetiDPCoarseProblemSolver(subdomains, reordering);

            public IFetiDPSubdomainMatrixManager CreateMatricesManager(ISubdomain subdomain)
                => new SkylineFetiDPSubdomainMatrixManager(subdomain, reordering);
        }
    }
}
