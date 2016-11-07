using System;
using System.Collections.Generic;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Solvers.Skyline
{
    public class SubdomainSkyline : ILinearSystem
    {
        private readonly Subdomain subdomain;
        private SkylineMatrix2D stiffnessMatrix;
        // REMOVE
        private SkylineMatrix2D stiffnessMatrixCopy;
        private Vector solution;

        public SubdomainSkyline(Subdomain subdomain)
        {
            this.subdomain = subdomain;
            solution = new Vector(subdomain.TotalDOFs);
        }
        
        #region ISolverSubdomain Members

        public int ID
        {
            get { return subdomain.ID; }
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
            get { return new Vector(subdomain.Forces); }
        }

        public IVector Solution
        {
            get { return solution; }
            set { solution = (Vector)value; }
        }

        //public void CloneMatrix()
        //{
        //    stiffnessMatrixCopy = (SkylineMatrix2D<double>)stiffnessMatrix.Clone();
        //}

        public IVector GetRHSFromSolution(IVector solution, IVector dSolution)
        {
            //// REMOVE
            //var forces = new double[solution.Length];
            //stiffnessMatrixCopy.Multiply(solution, forces);
            //return new Vector<double>(forces);

            //return subdomain.GetRHSFromSolution(solution, dSolution);
            throw new NotImplementedException("Check commented line of code above.");
        }

        public void SaveMaterialState()
        {
            subdomain.SaveMaterialState();
        }

        public void ClearMaterialStresses()
        {
            subdomain.ClearMaterialStresses();
        }

        public void SubdomainToGlobalVector(double[] vIn, double[] vOut)
        {
            foreach (int nodeID in subdomain.GlobalNodalDOFsDictionary.Keys)
            {
                Dictionary<DOFType, int> dofTypes = subdomain.NodalDOFsDictionary[nodeID];
                foreach (DOFType dofType in dofTypes.Keys)
                {
                    int localDOF = subdomain.NodalDOFsDictionary[nodeID][dofType];
                    int globalDOF = subdomain.GlobalNodalDOFsDictionary[nodeID][dofType];
                    if (localDOF > -1 && globalDOF > -1) vOut[globalDOF] += vIn[localDOF];
                }
            }
        }

        public void SubdomainToGlobalVectorMeanValue(double[] vIn, double[] vOut)
        {
            throw new NotImplementedException();
        }

        public void SplitGlobalVectorToSubdomain(double[] vIn, double[] vOut)
        {
            foreach (int nodeID in subdomain.GlobalNodalDOFsDictionary.Keys)
            {
                Dictionary<DOFType, int> dofTypes = subdomain.NodalDOFsDictionary[nodeID];
                foreach (DOFType dofType in dofTypes.Keys)
                {
                    int localDOF = subdomain.NodalDOFsDictionary[nodeID][dofType];
                    int globalDOF = subdomain.GlobalNodalDOFsDictionary[nodeID][dofType];
                    if (localDOF > -1 && globalDOF > -1) vOut[localDOF] = vIn[globalDOF];
                }
            }
        }

        #endregion

    }
}
