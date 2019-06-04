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
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: Add state checking for all the managed matrices. A state machine (using State pattern) should help.
namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1.Matrices
{
    public class DenseFeti1SubdomainMatrixManager : IFeti1SubdomainMatrixManager
    {
        private readonly SkylineAssembler assembler = new SkylineAssembler();
        private readonly SingleSubdomainSystem<SkylineMatrix> linearSystem;

        private SemidefiniteCholeskySkyline inverseKff;
        private Matrix inverseKii;
        private DiagonalMatrix inverseKiiDiagonal;
        private Matrix Kbb;
        private Matrix Kbi;

        public DenseFeti1SubdomainMatrixManager(ISubdomain subdomain)
        {
            this.linearSystem = new SingleSubdomainSystem<SkylineMatrix>(subdomain);
        }

        public ISingleSubdomainLinearSystem LinearSystem => linearSystem;

        public List<Vector> RigidBodyModes { get; private set; }

        public IMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement> elements, 
            IElementMatrixProvider matrixProvider)
        {
            return assembler.BuildGlobalMatrix(dofOrdering, elements, matrixProvider);
        }

        public (IMatrix Kff, IMatrixView Kfc, IMatrixView Kcf, IMatrixView Kcc) BuildGlobalSubmatrices(
            ISubdomainFreeDofOrdering freeDofOrdering, ISubdomainConstrainedDofOrdering constrainedDofOrdering, 
            IEnumerable<IElement> elements, IElementMatrixProvider matrixProvider)
        {
            return assembler.BuildGlobalSubmatrices(freeDofOrdering, constrainedDofOrdering, elements, matrixProvider);
        }

        public void Clear()
        {
            inverseKff = null;
            inverseKii = null;
            inverseKiiDiagonal = null;
            Kbb = null;
            Kbi = null;
            //linearSystem.Matrix = null; // DO NOT DO THAT!!! The analyzer manages that.
        }

        public void ExtractAndInvertKiiDiagonal(int[] internalDofs)
        {
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

        public void ExtractAndInvertKii(int[] internalDofs)
        {
            try
            {
                inverseKii = linearSystem.Matrix.GetSubmatrix(internalDofs, internalDofs);
                inverseKii.InvertInPlace();
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The free-free stiffness matrix of this subdomain has been overwritten and cannot be used anymore."
                    + "Try calling this method before factorizing/inverting it.");
            }
        }

        public void ExtractKbb(int[] boundaryDofs)
        {
            try
            {
                Kbb = linearSystem.Matrix.GetSubmatrix(boundaryDofs, boundaryDofs);
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The free-free stiffness matrix of this subdomain has been overwritten and cannot be used anymore."
                    + "Try calling this method before factorizing/inverting it.");
            }
        }

        public void ExtractKbiKib(int[] boundaryDofs, int[] internalDofs)
        {
            try
            {
                Kbi = linearSystem.Matrix.GetSubmatrix(boundaryDofs, internalDofs);
            }
            catch (MatrixDataOverwrittenException)
            {
                throw new InvalidOperationException(
                    "The free-free stiffness matrix of this subdomain has been overwritten and cannot be used anymore."
                    + "Try calling this method before factorizing/inverting it.");
            }
        }

        public void HandleDofOrderingWillBeModified() => assembler.HandleDofOrderingWillBeModified();

        public void InvertKff(double factorizationTolerance, bool inPlace)
        {
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
            return inverseKii * vector;
        }

        public Matrix MultiplyInverseKiiTimes(Matrix matrix)
        {
            if (inverseKii == null)
            {
                throw new InvalidOperationException("The inverse of the internal-internal stiffness submatrix"
                    + " 'inv(Kii)' of this subdomain must be calculated first.");
            }
            return inverseKii * matrix;
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
            if (Kbi == null)
            {
                throw new InvalidOperationException(
                    "The boundary-internal stiffness submatrix 'Kbi' of this subdomain must be calculated first.");
            }
            return Kbi * vector;
        }

        public Matrix MultiplyKbiTimes(Matrix matrix)
        {
            if (Kbi == null)
            {
                throw new InvalidOperationException(
                    "The boundary-internal stiffness submatrix 'Kbi' of this subdomain must be calculated first.");
            }
            return Kbi * matrix;
        }

        public Vector MultiplyKibTimes(Vector vector)
        {
            if (Kbi == null)
            {
                throw new InvalidOperationException(
                    "The boundary-internal stiffness submatrix 'Kbi' of this subdomain must be calculated first.");
            }
            return Kbi.Multiply(vector, true);
        }

        public Matrix MultiplyKibTimes(Matrix matrix)
        {
            if (Kbi == null)
            {
                throw new InvalidOperationException(
                    "The boundary-internal stiffness submatrix 'Kbi' of this subdomain must be calculated first.");
            }
            return Kbi.MultiplyRight(matrix, true);
        }

        public void SetSolutionVector(Vector solution) => linearSystem.SolutionConcrete = solution;

        public class Factory : IFeti1SubdomainMatrixManagerFactory
        {
            public IFeti1SubdomainMatrixManager CreateMatricesManager(ISubdomain subdomain)
                => new DenseFeti1SubdomainMatrixManager(subdomain);
        }
    }
}
