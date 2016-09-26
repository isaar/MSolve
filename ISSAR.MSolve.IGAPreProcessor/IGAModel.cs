using ISAAR.MSolve.Matrices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISSAR.MSolve.IGAPreProcessor
{
    public class IGAModel
    {
        private int numberOfDimensions;
        private int degreeKsi;
        private int degreeHeta;
        private int degreeZeta;
        private Vector<double> knotValueVectorKsi;
        private Vector<double> knotValueVectorHeta;
        private Vector<double> knotValueVectorZeta;
        private int totalDofs;
        private IList<ControlPoint> controlPoints = new List<ControlPoint>();
        private IList<Knot> knots = new List<Knot>();
        private IList<Element> elements = new List<Element>();


    }
}
