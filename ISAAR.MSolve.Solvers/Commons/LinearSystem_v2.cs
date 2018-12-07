using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: add state management
namespace ISAAR.MSolve.Solvers.Commons
{
    public abstract class LinearSystem_v2<TMatrix, TVector> : ILinearSystem_v2
        where TMatrix : IMatrixView //TODO: perhaps this should be IMatrix
        where TVector : IVector
    {
        protected LinearSystem_v2(ISubdomain_v2 subdomain)
        {
            this.Subdomain = subdomain;
        }

        public ISubdomain_v2 Subdomain { get; }

        public HashSet<ISystemMatrixObserver> MatrixObservers { get; } = new HashSet<ISystemMatrixObserver>();

        //TODO: I do not like this at all. However the providers are a mess right now, so this should be refactored after them.
        public bool IsMatrixOverwrittenBySolver { get; internal set; }

        IMatrixView ILinearSystem_v2.Matrix => Matrix;

        internal TMatrix Matrix { get; set; }

        IVector ILinearSystem_v2.RhsVector
        {
            get => RhsVector;
            set => RhsVector = (TVector)value;
        }

        internal TVector RhsVector { get; set; }

        IVectorView ILinearSystem_v2.Solution { get => Solution; }
        internal TVector Solution { get; set; }

        IVector ILinearSystem_v2.CreateZeroVector() => CreateZeroVector();

        public void SetMatrix(IMatrixView matrix)
        {
            foreach (var observer in MatrixObservers) observer.OnMatrixSetting();
            Matrix = (TMatrix)matrix;
            IsMatrixOverwrittenBySolver = false;
        }
        
        public abstract TVector CreateZeroVector();
        public abstract void GetRhsFromSubdomain();
    }
}
