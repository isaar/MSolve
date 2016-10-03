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
        private int id;
        private IList<Knot> knots;
        private IVector<int> connectivity;

        public NURBSElement3D(int id, IList<Knot> knots, IVector<int> connectivity)
        {
            this.id = id;
            this.knots = knots;
            this.connectivity = connectivity;
        }

        #region IISogeometricStructuralElement
        public ElementDimensions ElementDimensions
        {
            get
            {
               return ElementDimensions.ThreeD;
            }
        }

        public int ID
        {
            get
            {
                return this.id;
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
        #endregion
    }
}
