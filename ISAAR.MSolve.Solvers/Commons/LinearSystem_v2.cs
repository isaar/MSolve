using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: add state management
namespace ISAAR.MSolve.Solvers.Commons
{
    public class LinearSystem_v2<TMatrix, TVector>: ILinearSystem_v2
        where TMatrix:IMatrix 
        where TVector:IVector
    {
        public LinearSystem_v2(int id)
        {
            this.ID = id;
        }

        public int ID { get; }

        //TODO: this is error prone. This object should manage the state when clients read or modify the matrix.
        public bool IsMatrixFactorized { get; set; }
        public bool IsMatrixModified { get; set; } = true;

        IMatrix ILinearSystem_v2.Matrix
        {
            get => Matrix;
            set => Matrix = (TMatrix)value;
        }

        internal TMatrix Matrix { get; set; }

        IVector ILinearSystem_v2.RhsVector
        {
            get => RhsVector;
            set => RhsVector = (TVector)value; 
        }

        internal TVector RhsVector { get; set; }

        IVector ILinearSystem_v2.Solution { get => Solution; }
        internal TVector Solution { get; set; }
    }
}
