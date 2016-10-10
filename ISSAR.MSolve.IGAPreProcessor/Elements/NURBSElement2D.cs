using ISSAR.MSolve.IGAPreProcessor.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.PreProcessor.Elements.SupportiveClasses;
using ISAAR.MSolve.Matrices;

namespace ISSAR.MSolve.IGAPreProcessor.Elements
{
    public class NURBSElement2D:IIsogeometricStructuralElement
    {
        private int id;
        private IList<Knot> knots;
        private IVector<int> connectivity;
        private IGAModel model;
        private IList<GaussLegendrePoint3D> gaussPoints;
        
        public NURBSElement2D(int id,IGAModel model, IList<Knot> knots, IVector<int> connectivity)
        {
            this.id = id;
            this.knots = knots;
            this.connectivity = connectivity;
            this.model = model;

            
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

        private void CreateElementGaussPoints()
        {
            GaussLegendrePoint1D[] gaussPointsPerAxisKsi = 
                GaussQuadrature.GetGaussLegendrePoints(model.DegreeKsi);

            GaussLegendrePoint1D[] gaussPointsPerAxisHeta =
                GaussQuadrature.GetGaussLegendrePoints(model.DegreeHeta);

            IVector<double> coordinatesKsi = new Vector<double>(gaussPointsPerAxisKsi.Length);
            IVector<double> weightsKsi = new Vector<double>(gaussPointsPerAxisKsi.Length);
            for (int indexKsi = 0; indexKsi < gaussPointsPerAxisKsi.Length; indexKsi++)
            {
                coordinatesKsi[indexKsi] = 0.5 * (knots[0].Ksi + knots[2].Ksi+ (knots[2].Ksi - knots[0].Ksi) * gaussPointsPerAxisKsi[indexKsi].Coordinate);
                weightsKsi[indexKsi] = 0.5* ((knots[2].Ksi - knots[0].Ksi) * gaussPointsPerAxisKsi[indexKsi].WeightFactor);
            }

            IVector<double> coordinatesHeta = new Vector<double>(gaussPointsPerAxisHeta.Length);
            IVector<double> weightsHeta = new Vector<double>(gaussPointsPerAxisHeta.Length);
            for (int indexHeta = 0; indexHeta < gaussPointsPerAxisHeta.Length; indexHeta++)
            {
                coordinatesHeta[indexHeta] = 0.5 * (knots[0].Heta + knots[2].Heta + (knots[2].Heta - knots[0].Heta) * gaussPointsPerAxisHeta[indexHeta].Coordinate);
                weightsHeta[indexHeta] = 0.5 * ((knots[2].Heta - knots[0].Heta) * gaussPointsPerAxisHeta[indexHeta].WeightFactor);
            }



            



        }
    }
}
