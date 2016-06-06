using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.Solvers.Skyline
{
    public class SubdomainSkyline : ISolverSubdomain
    {
        private readonly Subdomain subdomain;
        private SkylineMatrix2D<double> stiffnessMatrix;
        // REMOVE
        private SkylineMatrix2D<double> stiffnessMatrixCopy;
        private Vector<double> solution;

        public SubdomainSkyline(Subdomain subdomain)
        {
            this.subdomain = subdomain;
            solution = new Vector<double>(subdomain.TotalDOFs);
        }
        
        #region ISolverSubdomain Members

        public int ID
        {
            get { return subdomain.ID; }
        }

        public IMatrix2D<double> Matrix
        {
            get { return stiffnessMatrix; }
            set 
            { 
                stiffnessMatrix = (SkylineMatrix2D<double>)value;
            }
        }

        public IVector<double> RHS
        {
            get { return new Vector<double>(subdomain.Forces); }
        }

        public IVector<double> Solution
        {
            get { return solution; }
            set { solution = (Vector<double>)value; }
        }

        //public void CloneMatrix()
        //{
        //    stiffnessMatrixCopy = (SkylineMatrix2D<double>)stiffnessMatrix.Clone();
        //}

        public IVector<double> GetRHSFromSolution(IVector<double> solution, IVector<double> dSolution)
        {
            //// REMOVE
            //var forces = new double[solution.Length];
            //stiffnessMatrixCopy.Multiply(solution, forces);
            //return new Vector<double>(forces);

            return subdomain.GetRHSFromSolution(solution, dSolution);
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
