using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.Matrices
{
    /// <summary>
    /// Dense format for Kii, Kbi/Kib, Kbb, Krr, Krc/Kcr, Kcc, KccStar and Skyline for Kff.
    /// Useful during prototyping and for debugging. For performance the other alternatives are probably better.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class DenseFetiDPSubdomainMatrixManager : IFetiDPSubdomainMatrixManager
    {
        private readonly SkylineAssembler assembler = new SkylineAssembler();
        private readonly SingleSubdomainSystem<SkylineMatrix> linearSystem;

        private DiagonalMatrix inverseKiiDiagonal;
        private Matrix inverseKii;
        private CholeskyFull inverseKrr;
        private Matrix Kbb;
        private Matrix Kbi;
        private Matrix Kcc;
        private Matrix KccStar; 
        private Matrix Krc;
        private Matrix Krr;

        public DenseFetiDPSubdomainMatrixManager(ISubdomain subdomain)
        {
            this.linearSystem = new SingleSubdomainSystem<SkylineMatrix>(subdomain);
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
        /// If the matrices stored in this object have already been calculated, they will be reused even if the original  
        /// free-free stiffness matrix has changed. To avoid that, this method must be called. 
        /// </summary>
        public void Clear()
        {
            inverseKii = null;
            inverseKiiDiagonal = null;
            inverseKrr = null;
            Kbb = null;
            Kbi = null;
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
            KccStar = Kcc - Krc.MultiplyRight(inverseKrr.SolveLinearSystems(Krc), true);
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
                inverseKii = Krr.GetSubmatrix(internalDofs, internalDofs);
                inverseKii.InvertInPlace();
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
                    //diagonal[i] = linearSystem.Matrix[idx, idx];
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
                Kbb = Krr.GetSubmatrix(boundaryDofs, boundaryDofs);
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
            if (Kbi != null) return;
            try
            {
                Kbi = Krr.GetSubmatrix(boundaryDofs, internalDofs);
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
            Krc = linearSystem.Matrix.GetSubmatrixFull(remainderDofs, cornerDofs);
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
            Krr = linearSystem.Matrix.GetSubmatrixFull(remainderDofs, remainderDofs);
        }

        public void HandleDofOrderingWillBeModified() => assembler.HandleDofOrderingWillBeModified();

        /// <summary>
        /// Will do nothing if it was already called. To perform this for a different stiffness matrix, first call 
        /// <see cref="Clear"/>.
        /// </summary>
        public void InvertKrr(bool inPlace)
        {
            if (inverseKrr != null) return;
            inverseKrr = Krr.FactorCholesky(inPlace);
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
            return inverseKii * vector;
        }

        public Matrix MultiplyInverseKiiTimes(Matrix matrix)
        {
            if (inverseKii == null)
            {
                throw new InvalidOperationException("The inverse of the internal_remainder - internal_remainder stiffness"
                    + " submatrix 'inv(Kii)' of this subdomain must be calculated first.");
            }
            return inverseKii * matrix;
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
            if (Kbi == null)
            {
                throw new InvalidOperationException("The boundary_remainder - internal_remainder stiffness submatrix"
                    + " 'Kbi' of this subdomain must be calculated first.");
            }
            return Kbi * vector;
        }

        public Matrix MultiplyKbiTimes(Matrix matrix)
        {
            if (Kbi == null)
            {
                throw new InvalidOperationException("The boundary_remainder - internal_remainder stiffness submatrix"
                    + " 'Kbi' of this subdomain must be calculated first.");
            }
            return Kbi * matrix;
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
            if (Kbi == null)
            {
                throw new InvalidOperationException("The internal_remainder - boundary_remainder stiffness submatrix"
                    + " 'Kib' of this subdomain must be calculated first.");
            }
            return Kbi.Multiply(vector, true);
        }

        public Matrix MultiplyKibTimes(Matrix matrix)
        {
            if (Kbi == null)
            {
                throw new InvalidOperationException("The internal_remainder - boundary_remainder stiffness submatrix"
                    + " 'Kib' of this subdomain must be calculated first.");
            }
            return Kbi.MultiplyRight(matrix, true);
        }

        public Vector MultiplyKrcTimes(Vector vector)
        {
            if (Krc == null)
            {
                throw new InvalidOperationException("The remainder-corner stiffness submatrix"
                    + " 'Krc' of this subdomain must be calculated first.");
            }
            return Krc * vector;
        }

        public void ReorderInternalDofs(FetiDPDofSeparator dofSeparator, ISubdomain subdomain)
        {
            // Do nothing, since the sparsity pattern is irrelevant for dense matrices.
        }

        public void ReorderRemainderDofs(FetiDPDofSeparator dofSeparator, ISubdomain subdomain)
        {
            // Do nothing, since the sparsity pattern is irrelevant for dense matrices.
        }

        public class Factory : IFetiDPSubdomainMatrixManagerFactory
        {
            public IFetiDPCoarseProblemSolver CreateCoarseProblemSolver(IReadOnlyList<ISubdomain> subdomains)
                => new DenseFetiDPCoarseProblemSolver(subdomains);

            public IFetiDPSubdomainMatrixManager CreateMatricesManager(ISubdomain subdomain)
                => new DenseFetiDPSubdomainMatrixManager(subdomain);
        }
    }
}
