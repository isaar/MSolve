using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.MaterialsTest
{
    public class ElementWithMaterial : IFiniteElement
    {
        private readonly IMaterialProperty materialProperty;
        private readonly IMaterialState[] materialStates;
        private readonly double[][] gaussPoints;

        public ElementWithMaterial(IMaterialProperty materialProperty)
        {
            this.materialProperty = materialProperty;
            // initialize gausspoints (we have coordinates in there)
            materialStates = new IMaterialState[gaussPoints.Length];
            for (int i = 0; i < gaussPoints.Length; i++)
                materialStates[i] = materialProperty.BuildMaterialState(gaussPoints[i]);
        }

        public double[] GetStiffnessMatrix()
        {
            //... integrate from materialStates
            materialStates[0].Stresses[0] = 0;
            return new double[0];
        }
    }
}
