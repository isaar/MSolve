using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: add state management
//TODO: should this and its implementations be internal? Analyzers and providers should work with ILinearSystem only. 
namespace ISAAR.MSolve.Solvers.Commons
{
    public abstract class LinearSystem_v2<TMatrix, TVector> : ILinearSystem_v2
        where TMatrix : class, IMatrix //TODO: perhaps this should be IMatrixView
        where TVector : class, IVector
    {
        private const int initialSize = int.MinValue;

        protected LinearSystem_v2(ISubdomain_v2 subdomain)
        {
            this.Subdomain = subdomain;
        }

        IMatrixView ILinearSystem_v2.Matrix => Matrix;

        public HashSet<ISystemMatrixObserver> MatrixObservers { get; } = new HashSet<ISystemMatrixObserver>();

        IVector ILinearSystem_v2.RhsVector
        {
            get => RhsVector;
            set
            {
                if (value.Length != this.Size)
                {
                    throw new NonMatchingDimensionsException("The provided vector does not match the dimensions or pattern of"
                        + " this linear system. Make sure that it is initialization was delegated to this linear system"
                        + " after the latest dof ordering.");
                }
                RhsVector = (TVector)value;
            }
        }

        public int Size { get; internal set; } = int.MinValue;

        public ISubdomain_v2 Subdomain { get; }

        IVectorView ILinearSystem_v2.Solution { get => Solution; }


        internal TMatrix Matrix { get; set; }

        internal TVector RhsVector { get; set; }

        internal TVector Solution { get; set; }

        public virtual void Clear()
        {
            foreach (var observer in MatrixObservers) observer.HandleMatrixWillBeSet();

            // Override the method if memory needs to be disposed in a more complicated way.
            RhsVector = null;
            Solution = null;
            Matrix = null;
        }

        IVector ILinearSystem_v2.CreateZeroVector()
        {
            if (Size == initialSize) throw new InvalidOperationException(
                "The linear system size must be set before creating vectors.");
            return CreateZeroVector();
        }

        public void SetMatrix(IMatrixView matrix)
        {
            if ((matrix.NumRows != this.Size) || (matrix.NumColumns != this.Size))
            {
                throw new NonMatchingDimensionsException("The provided matrix does not match the dimensions or pattern of"
                    + " this linear system. Make sure that it is initialization was delegated to this linear system"
                    + " after the latest dof ordering.");
            }
            foreach (var observer in MatrixObservers) observer.HandleMatrixWillBeSet();
            Matrix = (TMatrix)matrix;
        }

        internal abstract TVector CreateZeroVector();
    }
}
