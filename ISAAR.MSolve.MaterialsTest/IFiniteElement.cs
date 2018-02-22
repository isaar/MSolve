using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.MaterialsTest
{
    public interface IFiniteElement
    {
        //int ID { get; }
        //ElementDimensions ElementDimensions { get; }
        //IFiniteElementDOFEnumerator DOFEnumerator { get; set; }
        //IList<IList<DOFType>> GetElementDOFTypes(Element element);
        //bool MaterialModified { get; }
        //IMatrix2D<double> StiffnessMatrix(Element element);
        //IMatrix2D<double> MassMatrix(Element element);
        //IMatrix2D<double> DampingMatrix(Element element);
        //void ResetMaterialModified();
        //Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements);
        //double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements);
        //double[] CalculateForcesForLogging(Element element, double[] localDisplacements);
        //double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads);
        //void SaveMaterialState();
        //void ClearMaterialState();

        //void ClearMaterialStresses();
    }
}
