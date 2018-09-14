using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.Solvers.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.IGA.Writers
{
    public class VTUFileWriter
    {
        private readonly ILinearSystem LinearSystem;
        private readonly Model Model;
        private readonly int[] PointsPerElement;

        public VTUFileWriter(ILinearSystem linearSystem,Model model,int[] pointsPerElement)
        {
            LinearSystem = linearSystem;
            Model = model;
            PointsPerElement = pointsPerElement;
        }

        public void CreateTotalSolutionVector()
        {

        }


        private void CalculatePoints()
        {

        }

        private void CalculateDisplacements()
        {

        }

        public void CreateVTUDisplacementsFile()
        {

        }
    }
}
