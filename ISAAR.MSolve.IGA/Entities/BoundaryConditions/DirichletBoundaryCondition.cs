using System;
using System.Collections;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;

namespace ISAAR.MSolve.IGA.Entities.BoundaryConditions
{
    public class DirichletBoundaryCondition
    {
        private readonly Dictionary< int , Dictionary<IDofType,int>> dirichletDOFsDictionary = new Dictionary<int, Dictionary<IDofType, int>>();
        public double[] DirichletDisplacements { get; private set; }
        
        public DirichletBoundaryCondition(Model model, Func<double, double, double, double>dirichletValues, int[] dirichletSides)
        {
            //double[] loadVector = new double[model.ControlPointDOFsDictionary.Count];

            ArrayList dofs = new ArrayList();
            int numberOfValues = 0;
            foreach (int side in dirichletSides)
            {
                //TO DO: ask model for the dirichlet dofs for a specific side and add them to the dictionary
            }







        }
    }
}
