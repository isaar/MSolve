using ISSAR.MSolve.IGAPreProcessor.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISSAR.MSolve.IGAPreProcessor.Elements
{
    class NURBSElement3D : IIsogeometricStructuralElement
    {
        public ElementDimensions ElementDimensions
        {
            get
            {
                throw new NotImplementedException();
            }
        }

        public int ID
        {
            get
            {
                throw new NotImplementedException();
            }
        }

        public double[] CalculateForces(IGAElement element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(IGAElement element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public Tuple<double[], double[]> CalculateStresses(IGAElement element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public IMatrix2D<double> StiffnessMatrix(IGAElement element)
        {
            throw new NotImplementedException();
        }
    }
}
