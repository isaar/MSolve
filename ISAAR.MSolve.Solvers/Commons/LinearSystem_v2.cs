﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: add state management
namespace ISAAR.MSolve.Solvers.Commons
{
    public abstract class LinearSystem_v2<TMatrix, TVector> : ILinearSystem_v2
        where TMatrix : IMatrix
        where TVector : IVector
    {
        protected LinearSystem_v2(int id)
        {
            this.ID = id;
        }

        public int ID { get; }

        //TODO: this is error prone. This object should manage the state when clients read or modify the matrix.
        public bool IsMatrixFactorized { get; set; }
        public bool IsMatrixModified { get; set; } = true;
        public bool IsOrderModified { get; set; }

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

        IVectorView ILinearSystem_v2.Solution { get => Solution; }
        internal TVector Solution { get; set; }

        IVector ILinearSystem_v2.CreateZeroVector() => CreateZeroVector();

        ////TODO: this should probably be delegated to concrete classes, since it needs the current number of dofs, 
        ////      like CreateZeroVector(). Or it could be delegated to the solver.
        ////TODO: Checking the matrix columns may not be the best solution. In general, we need to check if the dofs have been  
        ////      changed (e.g.XFEM). If the dofs have been changed, but their number remains the same, we could still use the 
        ////      same memory.
        ////TODO: Take care of the case that this is called before a Solution vector has been assigned.
        //public void SetSolutionToZero() 
        //{ 
        //    // The order of the matrix might have been changed since the previous Solution vector had been created.
        //    if (Solution.Length == Matrix.NumColumns) Solution.Clear();
        //    else Solution = CreateZeroVector();
        //}

        public abstract TVector CreateZeroVector();
    }
}
