using System;
using System.Collections.Generic;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Skyline
{
    public class SkylineLinearSystem : ILinearSystem
    {
        private readonly int id;
        private readonly double[] forces;
        private SkylineMatrix2D stiffnessMatrix;
        // REMOVE
        private SkylineMatrix2D stiffnessMatrixCopy;
        private Vector solution;

        public SkylineLinearSystem(int id, double[] forces)
        {
            this.id = id;
            this.forces = forces;
            solution = new Vector(forces.Length);
        }

        #region ILinearSystem Members

        public int ID
        {
            get { return id; }
        }

        public IMatrix2D Matrix
        {
            get { return stiffnessMatrix; }
            set 
            { 
                stiffnessMatrix = (SkylineMatrix2D)value;
            }
        }

        public IVector RHS
        {
            get { return new Vector(forces); }
        }

        public IVector Solution
        {
            get { return solution; }
            set { solution = (Vector)value; }
        }

        #endregion

    }
}
