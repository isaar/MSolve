using ISSAR.MSolve.IGAPreProcessor.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices.Interfaces;

namespace ISSAR.MSolve.IGAPreProcessor.Elements
{
    public class NURBSElement2D:IIsogeometricStructuralElement
    {
        private int id;
        private IList<Knot> knots;
        private IVector<int> connectivity;
        
        public NURBSElement2D(int id, IList<Knot> knots, IVector<int> connectivity)
        {
            this.id = id;
            this.knots = knots;
            this.connectivity = connectivity;
        }
        #region IISogeometricStructuralElement
        public int ID
        {
            get
            {
                return this.id;
            }
        }

        public ElementDimensions ElementDimensions
        {
            get
            {
                return ElementDimensions.TwoD;
            }
        }

        public IMatrix2D<double> StiffnessMatrix(IGAElement element)
        {
            throw new NotImplementedException();
        }

        public Tuple<double[], double[]> CalculateStresses(IGAElement element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForces(IGAElement element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(IGAElement element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }
        
        #endregion
    }
}
