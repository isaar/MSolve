using ISAAR.MSolve.Matrices.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISSAR.MSolve.IGAPreProcessor.Interfaces
{
    public enum ElementDimensions
    {
        Unknown = 0,
        OneD = 1,
        TwoD = 2,
        ThreeD = 3
    }

    public interface IIsogeometricElement
    {
        int ID { get; }
        ElementDimensions ElementDimensions { get; }
       
        IMatrix2D<double> StiffnessMatrix(IGAElement element);
        Tuple<double[], double[]> CalculateStresses(IGAElement element, double[] localDisplacements, double[] localdDisplacements);
        double[] CalculateForces(IGAElement element, double[] localDisplacements, double[] localdDisplacements);
        double[] CalculateForcesForLogging(IGAElement element, double[] localDisplacements);
    }
}
